classdef Frame < handle
   properties
      empty
      model = []
      instance = []
      isEdge = false
      index = 0
      Preference = 0
      caseused = '-' %% Pred FMOd TbD Std
      interesting = false
      fittingScore = 0
      
      Area
      BoundingBox
      GTBoundingBox
      Centroid
      Start
      End
      Orientation
      PixelIdxList
      PixelList
      Boundary
      BoundaryXY
      Distances
      Radius
      Trajectory
      TrajectoryXY
      TrajectoryXYF
      TrajWeights
      Direction
      Length
      Speed
      SpeedFull
      f
      h
      hmask
      hi
      curve
      M
      T
      Ti
      coeff
      coeffi
      visim
      im_c
      bgr_c
      bb 
      mbb

      TP
      isFMO = false
      fmoCheck = false
      foundPrev = false
      foundNext = false
   end

   methods
      function this = Frame(region)
         if nargin == 0
            this.empty = true;
            return;
         end
         this.insert(region);
      end

      function [] = insert(this, region)
         for fn = fieldnames(region)'
            if isprop(this,fn{1})
               this.(fn{1}) = region.(fn{1});
            end
         end
         this.empty = false;
      end

      function [] = add(this, region)
         for fn = fieldnames(region)'
            if ~isempty(region.(fn{1}))
               this.(fn{1}) = region.(fn{1});
            end
         end
         this.empty = false;
      end

      function in = consist_color0(this, frames)
         in = Frame.consist_color(this, frames) < 0.3;
      end

      function in = consist_radius0(this, frames)
         in = Frame.consist_radius(this, frames) < 0.4;
      end

      function [FT, len] = remove_bg(this, im, bgr)
         M = APP.gen_ball(ceil(this.Radius));
         bb = this.BoundingBox;
         T = zeros(bb(3:4))';
         pxls = round(this.TrajectoryXY);
         this.BBTrajectory = [this.TrajectoryXY(1,:)-bb(1)+1; this.TrajectoryXY(2,:)-bb(2)+1];
         traj = sub2ind(size(T),pxls(2,:)-bb(2)+1,pxls(1,:)-bb(1)+1);
         T(traj) = 1;
         len = numel(traj);
         TM = conv2(T,M,'same')/len;
         FT = im.^(2.2) - (1-TM(:,:,[1 1 1])).*(bgr.^(2.2));
         % obj = sub2ind(size(T),this.PixelList(2,:)-bb(2)+1,this.PixelList(1,:)-bb(1)+1);
         % O = logical(zeros(size(TM))); O(obj) = 1;
         O = TM > 0;
         len = sum(O(:));
         FT(~O(:,:,[1 1 1])) = 0;
      end

      function add_crop(this, im1,im2)
         bbs = uint32(floor(this.BoundingBox));
         inc = 0.5;
         bbs(1:2) = bbs(1:2) - inc*bbs(3:4);
         bbs(3:4) = bbs(3:4) + 2*inc*bbs(3:4);

         bbs(1:2) = max(bbs(1:2), uint32([1 1]));
         bbs(3:4) = min([bbs(1:2)+bbs(3:4)], uint32([size(im2,2) size(im2,1)])) - bbs(1:2);

         this.CroppedBoundingBox = bbs;
         this.CroppedTrajectoryMask = logical(zeros(size(im2)));
         this.CroppedTrajectoryMask(this.Trajectory) = 1;
         this.CroppedTrajectoryMask = this.CroppedTrajectoryMask(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)));
         this.CroppedForeground = im2(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),:);
         this.CroppedBackground = im1(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),:);
      end

      function [] = add_info(this, sz)
         this.Trajectory = sub2ind(sz,this.TrajectoryXY(2,:),this.TrajectoryXY(1,:));
         this.PixelIdxList = sub2ind(sz,this.PixelList(2,:),this.PixelList(1,:));
         this.Length = numel(this.Trajectory);

         bin = logical(zeros(sz));
         bin(this.Trajectory) = 1;
         region = regionprops(bin,'Orientation');
         this.Orientation = region.Orientation;

         bin = logical(zeros(sz));
         bin(this.PixelIdxList) = 1;
         region = regionprops(bin,'BoundingBox');
         this.BoundingBox = region.BoundingBox;

         this.Edges = [this.First this.Last];
         this.Centroid = mean(this.PixelList')';
      end

      function [] = add_edges(this, frame0)
         if isempty(frame0.Last)
            [d,idx] = pdist2(this.Edges',frame0.Edges','euclidean','Smallest',1);
            [~,m] = min(d);
            frame0.Last = frame0.Edges(:,m);
            frame0.First = frame0.Edges(:,setdiff([1 2],m));
         else
            [d,idx] = pdist2(this.Edges',frame0.Last','euclidean','Smallest',1);
            [~,m] = min(d);
         end
         last = frame0.Last;
         this.First = this.Edges(:,idx(m));
         this.Last = this.Edges(:,setdiff([1 2],idx(m)));
         this.inter_traj = bwline2(this.First,last);
         this.inter_traj.pixels = this.inter_traj.pixels';
         this.inter_traj.extended = this.inter_traj.extended';
      end

      function [] = remove_edges(this)
         this.First = [];
         this.Last = [];
         this.inter_traj = [];
      end


      function [] = add_dist(this, sz)
         [y x] = ind2sub(sz,this.PixelIdxList);
         inl = x > 0 & y > 0 & x <= sz(2) & y <= sz(1);
         x = x(inl); y = y(inl); this.PixelIdxList = this.PixelIdxList(inl);
         this.PixelList = [x; y];
         this.add_boundary(sz);
         bin = logical(ones(sz));
         bin(this.PixelIdxList) = 0;
         dist_trans = bwdist(bin);
         dist_trans = dist_trans(this.PixelIdxList);
         this.Distances = dist_trans;
         this.Radius = max(this.Distances);

         normd = this.Distances/this.Radius;
         this.Trajectory = this.PixelIdxList(normd > 0.7);
         bin = logical(zeros(sz));
         bin(this.Trajectory) = 1;
         bin = bwmorph(bin,'thin',Inf);
         this.Trajectory = find(bin)';

         this.Length = numel(this.Trajectory);
         if isempty(this.Orientation)
            regions = regionprops(bin,'Orientation','Area');
            [~,m] = max([regions.Area]);
            this.Orientation = regions(m).Orientation;
         end
         this.Alpha = (2*this.Radius)/(this.Length + 2*this.Radius);

         [y x] = ind2sub(sz,this.Trajectory);
         this.TrajectoryXY = [x; y];

         [d,idx] = pdist2(this.TrajectoryXY',this.TrajectoryXY','euclidean','Largest',1);
         [~,m] = max(d);
         this.Area = numel(this.PixelIdxList);

         in = logical(zeros(sz));
         bin(this.PixelIdxList) = 1;
         regions = regionprops(bin,'BoundingBox','Centroid','Area');
         [~,maxarea] = max([regions.Area]);
         regions = regions(maxarea);
         this.BoundingBox = regions.BoundingBox;
         this.Centroid = regions.Centroid;
         if isempty(this.Direction)
            this.Direction = [1,2];
         end
         this.check1();
      end

      function [] = check1(this)
         assert(size(this.PixelIdxList,1) == 1,'Wrond dimension');
         assert(size(this.PixelList,1) == 2,'Wrond dimension');
         assert(size(this.Trajectory,1) == 1,'Wrond dimension');
         assert(size(this.TrajectoryXY,1) == 2,'Wrond dimension');
         assert(size(this.Boundary,1) == 1,'Wrond dimension');
         assert(size(this.BoundaryXY,1) == 2,'Wrond dimension');
         assert(size(this.Distances,1) == 1,'Wrond dimension');
      end

      function [] = check2(this)
         assert(all(size(this.Color) == [3 1] & size(this.MixedColor) == [3 1]),'Wrond dimension');
      end


      function [] = add_colors(this, back, front)
         offset = size(front,1)*size(front,2);
         ind = [this.Trajectory; ...
                this.Trajectory+offset; ...
                this.Trajectory+2*offset];
         colors = (front(ind) - (1-this.Alpha)*back(ind))/this.Alpha;
         
         if numel(ind) > 3
            this.MixedColor = mean(front(ind)')';
            this.Color = mean(colors')';
         else
            this.MixedColor = front(ind);
            this.Color = colors;
         end
         this.Color = max(min(this.Color, [1 1 1]'),[0 0 0]');
         this.check2();
      end

      function [] = save(this, file, im, frame0, dI, dI0, bin)
         offset = size(im,1)*size(im,2);
         bbs = this.BoundingBox;
         a = 30;
         bbs = bbs - [a a -2*a -2*a];
         bbs(1) = max(1,bbs(1)); bbs(2) = max(1,bbs(2));
         bbs(3) = min(size(im,2),bbs(1)+bbs(3)) - bbs(1);
         bbs(4) = min(size(im,1),bbs(2)+bbs(4)) - bbs(2);
         t = im;
         if nargin > 3 && ~frame0.empty
            a = 0;
            color0 = [0 0 255];
            bd0 = [frame0.Trajectory; frame0.Boundary];
            t(bd0) = a*t(bd0) + (1-a)*color0(1);
            t(bd0+offset) = a*t(bd0+offset) + (1-a)*color0(2);
            t(bd0+2*offset) = a*t(bd0+2*offset) + (1-a)*color0(3);
            % colort = [255 0 0];
            % traj = bwline2(frame0.Last,this.First);
            % bdt = [traj.pixels; traj.Boundary];
            % t(bdt) = a*t(bdt) + (1-a)*colort(1);
            % t(bdt+offset) = a*t(bdt+offset) + (1-a)*colort(2);
            % t(bdt+2*offset) = a*t(bdt+2*offset) + (1-a)*colort(3);
         end
         color = [0 255 0];
         bd = [this.Trajectory; this.Boundary];
         t(bd) = color(1);
         t(bd+offset) = color(2);
         t(bd+2*offset) = color(3);
         t = t(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),:);
         t = [t dI0(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),:); ...
               bin(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),[1 1 1]) ...
               dI(bbs(2):(bbs(2)+bbs(4)),bbs(1):(bbs(1)+bbs(3)),:)];
         imwrite(t,file);
      end

      function add_boundary(this,sz)
         x = this.PixelList(1,:);
         y = this.PixelList(2,:);
         bbs = [min(x) min(y) max(x) max(y)];
         bbs(3:4) = bbs(3:4) - bbs(1:2) + 1;
         img = logical(zeros(bbs(3),bbs(4)));
         t1 = sub2ind(bbs(3:4),x-bbs(1)+1,y-bbs(2)+1);
         img(t1) = 1;
         B = bwmorph(img','remove');
         [y x] = find(B);
         if size(B,1) == 1
            y = y'; x = x';
         end
         this.BoundaryXY = [x'+bbs(1)-1; y'+bbs(2)-1];
         this.Boundary = sub2ind(sz,this.BoundaryXY(2,:),this.BoundaryXY(1,:));
      end     

      function [] = show(this,color,lw,varargin)
         if nargin < 2
            color = 'r';
            if strcmp(this.caseused,'TbD') 
               color = 'g'; 
            elseif strcmp(this.caseused,'Pred')
               color = 'b';
            end
         end
         if nargin < 3
            lw = 10;
         end
         hold on;
         if this.TP == 2
            % rectangle('Position',this.bb); 
         end
         if ~isempty(this.Length)
            s = scatter(this.TrajectoryXY(1,:), ...
                  this.TrajectoryXY(2,:), ...
                  '.','LineWidth',lw,'MarkerEdgeColor',color,varargin{:});
            % if ~isempty(this.inter_traj) && ~isempty(this.inter_traj.pixels)
            %    s = scatter(this.inter_traj.pixels(1,:), ...
            %        this.inter_traj.pixels(2,:), ...
            %        '.','MarkerSize',lw,'MarkerEdgeColor',[0 0.5 1]);
            % end
            % s = scatter(this.BoundaryXY(1,:), ...
            %       this.BoundaryXY(2,:), ...
            %       '.','LineWidth',lw,'MarkerEdgeColor',color,varargin{:});
         end
         if ~isempty(this.End) && this.Preference == 2 && false
            dp = this.End - this.Start;
            st = this.Start;
            % st = this.Start + dp/2;
            % dp = dp / 2;
            quiver(st(1),st(2),dp(1),dp(2),0,'LineWidth',lw/2,'MaxHeadSize', 40, 'Color',color);
         end
      end

      function [] = show_number(this, num)
         if isempty(this.Centroid)
            this.Centroid = ((this.Start + this.End) / 2);
         end
         if isempty(this.Centroid)
            this.Centroid = this.bb(1:2) + (this.bb(3:4) / 2);
         end
         this.Centroid = double(this.Centroid);
         if ~isempty(this.Centroid)
            text(this.Centroid(1), this.Centroid(2), strcat('!',num2str(num),'!'), 'FontSize', 10);
         end
      end

      function [] = show_traj(this,color,lw,varargin)
         if nargin < 2
            color = 'b';
         end
         if nargin < 3
            lw = 10;
         end
         hold on;
         s = scatter(this.TrajectoryXY(1,:), ...
               this.TrajectoryXY(2,:), ...
               '.','LineWidth',lw,'MarkerEdgeColor',color,varargin{:});
      end

      function img = apply(this,img,color)
         if nargin < 3
            color = [0 0 1];
         end
         offset = size(img,1)*size(img,2);
         BW = logical(zeros(size(img,1),size(img,2)));
         BW(this.Boundary) = 1;
         BW(this.Trajectory) = 1;
         dist = bwdist(BW);
         BW2 = BW | (dist < 3);
         pxls = find(BW2);

         img(pxls) = color(1);
         img(pxls+offset) = color(2);
         img(pxls+2*offset) = color(3);
      end

   end

   methods(Static)
      function strct = forHonza(frames) 
         strct = [];
         for fr = [frames{:}]
            s.im = fr.im_c;
            s.bgr = fr.bgr_c;
            s.Ti = fr.Ti;
            s.T = fr.T;
            s.r = fr.Radius;
            strct = [strct s];
         end
      end

      function s = saveForHonza(im_c, bgr_c, T, r)
         global data;
         d.im = im_c;
         d.bgr = bgr_c;
         d.T = T;
         d.r = r;
         data = [data d];
      end

      function ori = check_ori(fmos, bin0, bin1)
         ori = 0;
         len = fmos.MajorAxisLength;
         ang = fmos.Orientation * pi / 180;
         pnts = fmos.PixelList;
         sz = size(bin0);
         d = [len:ceil(len/10):(2*len)];
         if isempty(d)
            return
         end
         cosine = cos(ang); sine = sin(ang);
         for i = 1:numel(d)
            pnts1 = round(bsxfun(@plus, pnts', d(i)*[cosine; -sine]));
            pnts0 = round(bsxfun(@plus, pnts', -d(i)*[cosine; -sine]));
            inl = pnts1(2,:) > 0 & pnts1(1,:) > 0 & pnts1(2,:) <= sz(1) & pnts1(1,:) <= sz(2);
            n_inlp(i) = sum(inl);
            pnts1 = sub2ind(sz, pnts1(2,inl), pnts1(1,inl));
            inl = pnts0(2,:) > 0 & pnts0(1,:) > 0 & pnts0(2,:) <= sz(1) & pnts0(1,:) <= sz(2);
            n_inlm(i) = sum(inl);
            pnts0 = sub2ind(sz, pnts0(2,inl), pnts0(1,inl));
            scorep(i) = min(sum(bin1(pnts1))/n_inlp(i), sum(bin0(pnts0))/n_inlm(i));
            scorem(i) = min(sum(bin1(pnts0))/n_inlm(i), sum(bin0(pnts1))/n_inlp(i));
         end
         min_th = Differential.min_cc/numel(fmos.PixelIdxList);
         scorep(scorep < min_th) = 0;
         scorem(scorem < min_th) = 0;
         p = max(scorep);
         m = max(scorem);
         if max(p,m) < 0.2
            return;
         elseif p > m
            ori = 1;
         else
            ori = -1;
         end
      end

      function d1 = consist_color(frame0, frames)
         [~,d] = knnMatch(frame0,frames,'mixcolor');
         [~,d0] = knnMatch(frame0,frames,'color');
         d1 = min([d d0]');
      end

      function d2 = consist_radius(frame0, frames)
         d2 = abs([frame0.Radius] - [frames.Radius])/frame0.Radius;
      end

      function is = is_connected(traj, sz)
         [y x] = ind2sub(sz,traj);
         bbs = [min(x) min(y) max(x) max(y)];
         bbs(3:4) = bbs(3:4) - bbs(1:2) + 1;
         img = logical(zeros(bbs(3),bbs(4)));
         t1 = sub2ind(bbs(3:4),x-bbs(1)+1,y-bbs(2)+1);
         img(t1) = 1;
         [~,num] = bwlabel(img);
         is = num == 1;
      end

      function frame = create(first, last, neighbours, model, sz)
         dx = last(1) - first(1);
         dy = last(2) - first(2);
         hyp = sqrt(dx^2 + dy^2);
         sine = dy/hyp;
         cosine = dx/hyp;
         region.TrajectoryXY = round(bsxfun(@plus, first, bsxfun(@mtimes, [1:hyp],[cosine; sine])));
         d = pdist2(region.TrajectoryXY',neighbours','euclidean','Smallest',1);
         inl = round(neighbours(:,d < model.Radius));
         dist = d(d < model.Radius);
         dist = dist(inl(2,:) > 0 & inl(1,:) > 0 & inl(2,:) <= sz(1) & inl(1,:) <= sz(2));
         inl = inl(:,inl(2,:) > 0 & inl(1,:) > 0 & inl(2,:) <= sz(1) & inl(1,:) <= sz(2));
         inrange = region.TrajectoryXY(2,:) > 0 & region.TrajectoryXY(1,:) > 0 & region.TrajectoryXY(2,:) <= sz(1) & region.TrajectoryXY(1,:) <= sz(2);
         region.TrajectoryXY = region.TrajectoryXY(:,inrange);
         region.PixelList = inl;
         region.Distances = model.Radius - dist;
         region.First = first;
         region.Last = last;
         region.Radius = model.Radius;

         if ~any(inrange)
            frame = Frame();
            return;
         elseif ~all(inrange)
            [d,idx] = pdist2(region.TrajectoryXY',region.TrajectoryXY','euclidean','Largest',1);
            [~,m] = max(d);
            region.Edges = region.TrajectoryXY(:,[m idx(m)]);
            d2first = sum((region.Edges - [first first]).^2);
            [~,ind] = sort(d2first);
            region.First = region.Edges(:,ind(1));
            region.Last = region.Edges(:,ind(2));
         end
         region.PixelIdxList = sub2ind(sz,region.PixelList(2,:),region.PixelList(1,:));
         frame = Frame(region);
         frame.add_dist(sz);
      end

      function bbs = combine_bbs(bbs1,bbs2)
         if isempty(bbs1)
            bbs = bbs2;
            return;
         elseif isempty(bbs2)
            bbs = bbs1;
            return
         end
         bbs(1) = min(bbs1(1),bbs2(1));
         bbs(2) = min(bbs1(2),bbs2(2));
         bbs(3) = max(bbs1(1)+bbs1(3),bbs2(1)+bbs2(3)) - bbs(1);
         bbs(4) = max(bbs1(2)+bbs1(4),bbs2(2)+bbs2(4)) - bbs(2);
      end
   end

end
