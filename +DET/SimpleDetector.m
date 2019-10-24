classdef SimpleDetector < handle
%% Detector with the following assumptions:
%% 1) There is exactly one FMO in each frame
%% 2) There is no other moving object in the whole sequence
%% 3) Background is given, the camera is static
%% 4) Trajectory within one frame in linear -> bounces not detected
	properties
		BGR = []
		noise = 20/255
		linear_fit = true;
		linear_fit_init = true;
		do_deblur = true;
	end

	methods
		function this = SimpleDetector(varargin)
			this = cmp_argparse(this, varargin);
		end

		function fr = detect(this, IM)
			params_psffit = []; params_psffit.seqransac_lin_max_samples = 7140; params_psffit.ransac_quad_max_samples = 4000;

			delta = sum(abs(IM - this.BGR),3) > this.noise;
			delta = imfill(delta,'holes');
			regions = regionprops(delta, 'Area', 'PixelIdxList');
			if isempty(regions), error('Empty regions!'); end
			[~,ir] = max([regions.Area]); ir = ir(1);
			rgn = regions(ir);
			BW = zeros(size(delta));
			BW(rgn.PixelIdxList) = 1;
			CH = bwconvhull(BW);
			rgn = regionprops(CH, 'PixelList', 'Area', 'PixelIdxList', 'BoundingBox', 'Centroid');
			if numel(rgn) ~= 1, error('Not 1 region'); end

			dt = round(bwdist(~delta)); 
			lm = logical(dt >= imdilate(dt, [1 1 1; 1 0 1; 1 1 1])) & (dt > 1.5); 

			rgn.Trajectory = rgn.PixelIdxList(lm(rgn.PixelIdxList));
			Distances = dt(rgn.PixelIdxList);
			DistancesTrajectory = dt(rgn.Trajectory);
			rgn.Radius = max(DistancesTrajectory);
			maxs = DistancesTrajectory - 0.95*rgn.Radius >= -1;
			rgn.Trajectory = rgn.Trajectory(maxs);
			DistancesTrajectory = dt(rgn.Trajectory);

			[y,x] = ind2sub(size(IM), rgn.Trajectory);
			rgn.TrajectoryXY = [y x]';
			
			rgn.bb = bbextend(rgn.BoundingBox, 1.5*rgn.Radius + 5, size(IM));

			xy = rgn.TrajectoryXY - rgn.bb(2:-1:1)' + 1;
			
			rgn.hi = zeros(fliplr(rgn.bb(3:4) - rgn.bb(1:2)+1));
			rgn.hi(sub2ind(size(rgn.hi), xy(1,:), xy(2,:))) = DistancesTrajectory;

			if this.linear_fit_init
				try
					[rgn.Ti, rgn.coeffi] = psffit_linear(rgn.hi, 0.2);
				catch
					[rgn.coeffi, rgn.Ti] = psffit(rgn.hi, params_psffit);
				end
			else
				[rgn.coeffi, rgn.Ti] = psffit(rgn.hi, params_psffit);
			end

			rgn.h = rgn.hi;
			rgn.T = rgn.Ti;
			rgn.coeff = rgn.coeffi;

			M = diskMask([], rgn.Radius+7);
			im_c = IM(rgn.bb(2):rgn.bb(4),rgn.bb(1):rgn.bb(3),:);
			bgr_c = this.BGR(rgn.bb(2):rgn.bb(4),rgn.bb(1):rgn.bb(3),:);
			
			if this.do_deblur
				try
					[rgn.h rgn.f rgn.m] = blindLoop_FMH_TbD(im_c, bgr_c, repmat(M,[1 1 3]), M, [], [], 'lambda_R', 1e-2,'maxiter',300);
					bbox = CH(rgn.bb(2):rgn.bb(4),rgn.bb(1):rgn.bb(3),:);
					rgn.h(~bbox) = 0;
					if this.linear_fit
						[rgn.T, rgn.coeff] = psffit_linear(rgn.h, 0.2);
					else
						[rgn.coeff,rgn.T,val] = psffit(rgn.h, params_psffit);
					end
				catch
					disp('Error in blindloop')
				end
			end

			[y,x] = find(rgn.T > 0);
			rgn.TrajectoryXYF = [y'; x'] + rgn.bb(2:-1:1)' - 1;
			rgn.TrajectoryXY = rgn.TrajectoryXYF(2:-1:1,:);
			rgn.Length = numel(rgn.Trajectory);

			fr = VID.Frame(rgn);
			fr.im_c = im_c;
			fr.bgr_c = bgr_c;
			fr.caseused='TbD';
			fr.fittingScore = [];

			% if this.frmid > 20, keyboard; end

		end

	end

		
end