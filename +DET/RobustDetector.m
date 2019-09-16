classdef RobustDetector < handle
	properties
		BGR = [] % if static background
		B = []
		im0 = []
		im00 = []
		im000 = []
		noise_t = [0.06 0.06 0.06]
		noise = 5/255
		maxAllowedExp = 1.5
		Size = []
		Size_r = []
		scale_factor = 1
		counter = 0
		frame0 = []
		gm = 2.2
		score_th = 0.4
		corr_th = 0.95
		hs_th = 0.6
		min_f_th = -0.1
		min_radius = 3

		num_appear = 2
		do_tbd = true
		do_em = true
		do_stabilize = false
		linear_fit = false
		get_visim = false

		fast_version = false
	end

	methods
		function this = RobustDetector(sz, varargin)
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

		function frameOut = detect(this, IM)
			this.counter = this.counter + 1;
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
				if regions(k).Area < 300, continue; end
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
				regions(k).bb = bbextend(regions(k).BoundingBox, 0.8*regions(k).Radius, this.Size_r);

				regions(k).Length = numel(regions(k).Trajectory);
				if regions(k).Length <= 5, continue; end
				regions(k).DistancesTrajectory = dt(regions(k).Trajectory);

				[xs,ys] = ind2sub(this.Size_r, regions(k).Trajectory);
				regions(k).TrajectoryXY = [xs ys]';
				xy = regions(k).TrajectoryXY - regions(k).bb(2:-1:1)' + 1;
				
				% w = regions(k).DistancesTrajectory' / max(regions(k).DistancesTrajectory);
				% [regions(k).TrajectoryXYF, regions(k).TrajWeights] = fittraj(regions(k).TrajectoryXY,w,1);
				regions(k).Ti = zeros(fliplr(regions(k).bb(3:4) - regions(k).bb(1:2)+1));
				regions(k).Ti(sub2ind(size(regions(k).Ti), xy(1,:), xy(2,:))) = regions(k).DistancesTrajectory;
				% try
					if this.linear_fit
						[regions(k).Ti, regions(k).coeff] = psffit_linear(regions(k).Ti, 0.2);
					else
						[regions(k).Ti, regions(k).coeff] = psffit2(regions(k).Ti, 0.2);
					end
				% catch
				% 	disp('Error in fitting. [RobustDetector]');
				% 	continue;
				% end
				% if regions(k).Area  > 1500, keyboard; end

				regions(k).coeffi = regions(k).coeff;
				[y,x] = find(regions(k).Ti > 0);
				regions(k).TrajectoryXYF = [y'; x'] + regions(k).bb(2:-1:1)' - 1;

				if any(isnan(regions(k).TrajectoryXYF(:))), continue; end
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
				minRmse = Inf;
				indb = -1;

				if this.num_appear > 1
					for kfr0 = 1:numel(this.frame0)
						fr0 = this.frame0(kfr0);
						% if ~fr0.isFMO && ~fr0.isEdge, continue; end

					  	d = norm(fr.Centroid - fr0.Centroid);
		                if(d > this.maxAllowedExp*(fr.Length+2*fr.Radius)), continue; end
		                if(d < 0.5*(fr.Length+2*fr.Radius)), continue; end
		                if abs(fr.Length - fr0.Length)/min([fr0.Length fr.Length]) > 3, continue; end
		                if abs(fr.Radius - fr0.Radius)/min([fr0.Radius fr.Radius]) > 0.8, continue, end
		                xxxs = [fr0.TrajectoryXYF(1,:) fr.TrajectoryXYF(1,:)]';
		                yyys = [fr0.TrajectoryXYF(2,:) fr.TrajectoryXYF(2,:)]';
		                [~,g1] = fit(yyys,xxxs,'poly2');
						[~,g2] = fit(xxxs,yyys,'poly2');

		                rmse = min([g1.rmse g2.rmse]);
		                if rmse < minRmse 
		                	minRmse = rmse;
		                	indb = kfr0;
		                end
					end
					if minRmse < 1.5*fr.Radius 
						fr.foundPrev = true;
						fr.model = [fr.model indb];
					end
				end

				%%%%% DEBLURRING %%%%%%%%%%%%%%%
				doDeblur = (~fr.isEdge && fr.foundPrev) || this.num_appear == 1;
				if doDeblur
					sc = this.show_deblur(fr, im, bgr);
					if sc >= this.score_th 
						for frim = fr.model
							this.frame0(frim).foundNext = true;
						end
						fr.isFMO = true;
					end
				end
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			end
			
			%%%%% TbD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if this.do_tbd
				for fr0 = this.frame0
					if ~(fr0.isFMO && ~fr0.foundNext && fr0.foundPrev), continue; end
					[frmNext] = predict(this, fr0, im, bgr);
					global debugmode; 
					if ~isempty(debugmode), keyboard; end
					if ~isempty(frmNext)
						frame = [frame frmNext];
						fr0.foundNext = true;
					end
				end
			end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			frameOut = Frame.empty;
			if this.num_appear > 1
				if ~isempty(this.frame0)
					frameOut = this.getFrames(im, bgr);
					global pdfim;
					for fr = frameOut
						pdfim{end+1}.im = fr.visim;
						pdfim{end}.errbgr = 0; pdfim{end}.errim = 0;
					end
				end
			elseif this.num_appear == 1
				frameOut = frame([frame.isFMO]);
			end
			% if ~isempty(frame), keyboard; end

			this.frame0 = frame;
			this.proceed(im);
		end

		function [frameOut] = getFrames(this, im, bgr)
			frameOut = Frame.empty;
			if isempty(this.frame0), return; end
			conf = [this.frame0.foundPrev] | [this.frame0.foundNext];
			frameOut = this.frame0(conf);
			fmos = [frameOut.isFMO];
			for fri = 1:numel(fmos)
				if fmos(fri), continue; end
				sc = this.show_deblur(frameOut(fri), this.im0, bgr);
				if sc >= this.score_th 
					fr.isFMO = true;
					fmos(fri) = true;
				end		
			end
			frameOut = frameOut(fmos);
		end

		function [frmNext] = predict(this, frm, im, bgr)
			frmNext = Frame.empty;
			ex = max(frm.BoundingBox(3:4));
			bb = round(bbextend(frm.BoundingBox, ex, this.Size_r));
			im_c = im(bb(2):bb(4),bb(1):bb(3),:);
			bgr_c = bgr(bb(2):bb(4),bb(1):bb(3),:);
			binim = zeros(size(im,1),size(im,2));
			binim(frm.PixelIdxList) = 1;
			binim = binim(bb(2):bb(4),bb(1):bb(3));

			bgr_c_nogm = bgr_c.^(1/this.gm);
			im_c_nogm = im_c.^(1/this.gm);

			[R rsize] = createRotationMatrices(size(frm.M), 0);
			h = estimateH_rot(im_c_nogm, bgr_c_nogm, frm.f, frm.M, [], [], R, rsize, [], 'alpha', 1, 'maxiter', 40, 'angles',0);
			% h = estimateH_motion(im_c_nogm, bgr_c_nogm, frm.f, frm.M, [], [], [], 'alpha', 1, 'maxiter', 40);
			h(h<0) = 0;

			if sum(h(:)) < 1-this.hs_th || sum((h(find(binim)))) > 1+this.hs_th 
				return; 
			end
			try
				[T, coeffi] = psffit2(h, [0.3 0.1]);
			catch
				disp('Error in fitting. [DET.RobustDetector]');
				return;
			end
			pr = this.probH(h,T);
			if pr < this.score_th || pr > 2-this.score_th
				return; 
			end

			% f = estimateF_clr(im_c_nogm, bgr_c_nogm, T, frm.M, [], [], 'alpha', 2^-8, 'lp', 2, 'maxiter', 10, 'f_max', 10);
			f = estimateF_rot(im_c_nogm, bgr_c_nogm, T, frm.M, [], R, rsize, [], 'alpha', 2^-8, 'maxiter', 10, 'f_max', 10);
			rgn.Radius = frm.Radius;
			[y,x] = find(T > 0);
			rgn.TrajectoryXYF = [y'; x'] + bb(2:-1:1)' - 1;
			TM = conv2(T,frm.M,'same');
			[xx yy] = find(TM);
			rgn.PixelList = [yy'; xx'] + bb(1:2)' - 1;
			rgn.PixelIdxList = sub2ind(size(im),rgn.PixelList(2,:),rgn.PixelList(1,:));;

			rgn.TrajectoryXY = rgn.TrajectoryXYF(2:-1:1,:);
			rgn.f = f; 
			rgn.T = T; 
			rgn.Length = size(rgn.TrajectoryXYF,2)/2; 
			rgn.M = frm.M;
			rgn.Centroid = fliplr(mean(rgn.TrajectoryXYF'));
			rgn.coeffi = coeffi;
			rgn.coeff = coeffi;
			rgn.bb = bb;
			
		 	d = norm(frm.Centroid - rgn.Centroid);
            if( d > this.maxAllowedExp*(frm.Length+2*frm.Radius) ...
            	|| d < 0.3*(frm.Length) ...
            	|| abs(rgn.Length - frm.Length)/min([frm.Length rgn.Length]) > 0.5)
            	return; 
            end
	                
			rgn.BoundingBox = traj2bb(rgn.TrajectoryXYF, rgn.Radius);
			frmNext = Frame(rgn); frmNext.isFMO = true; frmNext.foundPrev = true;
			frmNext.im_c = im_c; 
			frmNext.bgr_c = bgr_c; 
			frmNext.fmoCheck = true;
			frmNext.visim = [];
			if this.get_visim
				frmNext.visim = genvisim(rgn.T,rgn.T,rgn.f,rgn.M,h,bgr_c,im_c,zeros(size(im_c)),this.gm);
			end
		end


		function [score] = show_deblur(this, fr, im, bgr)
			score = [];
			bb = fr.bb;
			im_c = im(bb(2):bb(4),bb(1):bb(3),:);
			bgr_c = bgr(bb(2):bb(4),bb(1):bb(3),:);
			binim = zeros(size(im,1),size(im,2));
			binim(fr.PixelIdxList) = 1;
			binim = binim(bb(2):bb(4),bb(1):bb(3));
			binim3 = repmat(binim,[1 1 3]);
			
			% Ti = p2T(fr.TrajectoryXYF, fr.TrajWeights, bb);

			fr.M = double(diskMask([],ceil(fr.Radius+1)));
			
			if this.do_em && sum(binim(:))/numel(im(:,:,1)) < 0.08
				% if this.fast_version
				% 	[h,fr.T,fr.f,fr.coeff,score,fr.fittingScore] = this.fitEMfast(im_c, bgr_c, fr.M, fr.Ti, binim);
				% else
				[h,fr.T,fr.f,fr.coeff,score,fr.fittingScore] = this.fitEM(im_c, bgr_c, fr.M, fr.Ti, binim);
				% end
			else
				h = fr.Ti; fr.T = fr.Ti;  fr.f = double(repmat(fr.M, [1 1 3]));
				try
					[~, fr.coeff] = psffit2(h, [0.3 0.1]);
				catch
					score = 0;
					return;
				end
				score = 1;  fr.fittingScore = 1;
			end
			score = min([score this.probT(fr.Ti, fr.T)]);

			if score > this.score_th && this.probT(fr.T,fr.Ti) > this.score_th
				[y,x] = find(fr.T > 0);
				fr.TrajectoryXY = [x'; y'] + bb(1:2)' - 1;
			end
			fr.visim = [];
			if this.get_visim
				fr.visim = genvisim(fr.Ti,fr.T,fr.f,fr.M,h,bgr_c,im_c,binim3,this.gm);
			end
			fr.im_c = im_c; fr.bgr_c = bgr_c; 
			fr.fmoCheck = true;
			fr.h = h;
			global debugmode; 
			if ~isempty(debugmode) && score >= 0, keyboard; end
		end

		function [h,T,f,coeff,prT0,fittingScore] = fitEMfast(this, im_c, bgr_c, M, T, binim)
			bgr_c_nogm = bgr_c.^(1/this.gm);
			im_c_nogm = im_c.^(1/this.gm);
			prT0 = 0;
			fittingScore = 0;
			coeff = [];

			T = T / sum(T(:));
			[h,f,m,~] = blindLoop_FMH_TbD(im_c_nogm, bgr_c_nogm, [], M, T, [], [], 'maxiter', 15,'rel_tol', 1e-2, 'cg_maxiter', 5);
			f(f<0) = 0; 

			% if sum((h(:))) < 1-this.hs_th, return; end
			if sum((h(find(binim)))) > 1+this.hs_th, return; end
			if this.probH(h,T) < 0.3, return; end
			[T0, coeff, fittingScore] = FIT.fitting(h,0.5, 2);
			prT0 = this.probT(T0, T);
			if prT0 < this.score_th, return; end
			T = T0;
		end

		function [h,T,f,coeff,prT0,fittingScore] = fitEM(this, im_c, bgr_c, M, T, binim)
			[R rsize] = createRotationMatrices(size(M), 0);
			bgr_c_nogm = bgr_c.^(1/this.gm);
			im_c_nogm = im_c.^(1/this.gm);
			prT0 = 0;
			fittingScore = 0;
			coeff = [];
			for ems = 1:3
				h = T;
				try
					f = estimateF_rot(im_c_nogm, bgr_c_nogm, h, M, [], R, rsize, [], 'alpha', 2^-8, 'maxiter', 10, 'f_max', 10);
				catch
					f = repmat(M, [1 1 3]);
					break;
				end
				f(f<0) = 0;
				% h = h / sum(h(:));
				h = estimateH_rot(im_c_nogm, bgr_c_nogm, f, M, [], [], R, rsize, [], 'alpha', 1, 'maxiter', 40, 'angles', 0);
				% h = estimateH_motion(im_c_nogm, bgr_c_nogm, f, M, [], [], [], 'alpha', 1, 'maxiter', 40);
				h(h<0) = 0;
				if sum((h(:))) < 1-this.hs_th, break; end
				if sum((h(find(binim)))) > 1+this.hs_th, break; end
				if this.probH(h,T) < 0.3, break; end
				% [T0, coeff] = psffit2(h, [0.3 0.1]); 
				[T0, coeff, fittingScore] = FIT.fitting(h,0.5, 2);
				prT0 = this.probT(T0, T);
				if prT0 < this.score_th
					break
				end
				T = T0;
				if prT0 > this.corr_th
					break
				end
			end

			if prT0 >= this.score_th
				try
					f = estimateF_rot(im_c_nogm, bgr_c_nogm, h, M, [], R, rsize, [], 'alpha', 2^-8, 'maxiter', 10, 'f_max', 10);
				catch
					prT0 = 0;
				end
			end
		end

		function prob = probT(this, T0, T) % pr(T0 | T)
			T0s = imdilate(T0>0,strel('square',3));
			prob = sum(T(T0s));
		end

		function prob = probH(this, h, T) % pr(H | T)
			T0s = imdilate(T>0,strel('square',3));
			prob = sum(h(T0s));
			% prob = sum(h(T0s))/sum(h(find(binim)));
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