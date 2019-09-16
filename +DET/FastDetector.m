classdef FastDetector < handle
	properties
		BGR = [] % if static background
		B = []
		im0 = []
		im00 = []
		im000 = []
		noise = 5/255
		max_allowed_exp = 1.5
		max_size = 300
		Size = []
		Size_r = []
		scale_factor = 1
		score_th = 0.4
		corr_th = 0.95
		hs_th = 0.6
		min_f_th = -0.1
		min_radius = 3
		do_stabilize = false
	end

	methods
		function this = FastDetector(sz, varargin)
			this.Size = sz;
			this.Size_r = round(this.Size * this.scale_factor);
			this = cmp_argparse(this, varargin);
		end

		function [] = nodetect(this, IM)
			im = imresize(IM, this.scale_factor);
			bgr = this.get_bgr(im);
			this.frame0 = Frame.empty;
			this.proceed(im);
		end

		function frame = detect(this, IM)
			im = imresize(IM, this.scale_factor);
			bgr = this.get_bgr(im);
			if this.do_stabilize
            	bgr = stabilize(im, bgr);
            end
			delta = this.binarize(abs(im - bgr));
			regions = regionprops(delta, 'PixelList', 'Area', 'PixelIdxList', 'BoundingBox', 'Centroid');
			delta_aug = [zeros(1,size(delta,2)); delta; zeros(1,size(delta,2))];
			delta_aug = [zeros(size(delta_aug,1),1) delta_aug zeros(size(delta_aug,1),1)];
			regions = regions([regions.Area] > 10);
			dt = round(bwdist(~delta_aug)); dt = dt(2:(end-1),2:(end-1),:);
			lm = logical(dt >= imdilate(dt, [1 1 1; 1 0 1; 1 1 1])) & (dt > 1.5); 

			frame = Frame.empty;
			for k = 1:numel(regions)
				if regions(k).Area < this.max_size, continue; end
				cnt = this.countEdge(regions(k).PixelList);
				if cnt > max(regions(k).BoundingBox(3:4))/2, continue; end
				regions(k).BoundingBox = round(regions(k).BoundingBox);
				
				regions(k).Trajectory = regions(k).PixelIdxList(lm(regions(k).PixelIdxList));
				regions(k).Distances = dt(regions(k).PixelIdxList);
				regions(k).DistancesTrajectory = dt(regions(k).Trajectory);
				regions(k).Radius = prctile(regions(k).DistancesTrajectory, 90);
				if(regions(k).Radius < this.min_radius), continue; end
				maxs = regions(k).DistancesTrajectory - 0.9*regions(k).Radius >= -1;
				regions(k).Trajectory = regions(k).Trajectory(maxs);
				regions(k).bb = bbextend(regions(k).BoundingBox, 1.5*regions(k).Radius + 5, this.Size_r);

				regions(k).Length = numel(regions(k).Trajectory);
				if regions(k).Length <= 5, continue; end
				regions(k).DistancesTrajectory = dt(regions(k).Trajectory);

				[xs,ys] = ind2sub(this.Size_r, regions(k).Trajectory);
				regions(k).TrajectoryXY = [xs ys]';
				xy = regions(k).TrajectoryXY - regions(k).bb(2:-1:1)' + 1;
				
				regions(k).Ti = zeros(fliplr(regions(k).bb(3:4) - regions(k).bb(1:2)+1));
				regions(k).Ti(sub2ind(size(regions(k).Ti), xy(1,:), xy(2,:))) = regions(k).DistancesTrajectory;
				[regions(k).Ti, regions(k).coeff] = psffit_linear(regions(k).Ti, 0.2);

				regions(k).coeffi = regions(k).coeff;
				[y,x] = find(regions(k).Ti > 0);
				regions(k).TrajectoryXYF = [y'; x'] + regions(k).bb(2:-1:1)' - 1;

				if any(isnan(regions(k).TrajectoryXYF)), continue; end
				if isempty(regions(k).TrajectoryXYF), continue; end
				try
					inds = sub2ind(this.Size_r, regions(k).TrajectoryXYF(1,:),regions(k).TrajectoryXYF(2,:));
				catch, continue; end
				if (sum(dt(inds) > 0) / numel(inds)) < 1, continue; end

				regions(k).TrajectoryXY = regions(k).TrajectoryXYF(2:-1:1,:);
				regions(k).hi = regions(k).Ti;

				frame = [frame regions(k)];
				fr = frame(end);
				fr.isEdge = (cnt > 0);
				fr.im_c = im(fr.bb(2):fr.bb(4),fr.bb(1):fr.bb(3),:);
				fr.bgr_c = bgr(fr.bb(2):fr.bb(4),fr.bb(1):fr.bb(3),:);
				fr.caseused='FMOd';
				fr.fittingScore = [];
			end
			this.proceed(im);
		end


		function bin = binarize(this, delta)
			bin = sum(delta,3) > this.noise;
		end

		function [] = proceed(this, im)
			this.im000 = this.im00; 
			this.im00 = this.im0;
			this.im0 = im;
		end

		function cnt = countEdge(this, pxls)
			cnt = sum(pxls(:,1) == this.Size(2)) + sum(pxls(:,2) == this.Size(1));
			cnt = cnt + sum([pxls(:)] == 1);
		end

		function bgr = get_bgr(this, im)
			if ~isempty(this.BGR)
				bgr = this.BGR;
			elseif isempty(this.im0) 
				bgr = im;
			elseif isempty(this.im00)
				bgr = this.im0;
			elseif isempty(this.im000)
				bgr = fast_median(im,this.im0,this.im00);
			else
				% bgr = fast_median(this.im0,this.im00,this.im000);
				bgr = median(cat(4,this.B,this.im0,this.im00,this.im000,im),4);
			end
			this.B = bgr;
		end
	end

		
end