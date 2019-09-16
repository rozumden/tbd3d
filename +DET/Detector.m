classdef Detector < DET.IDetector
	properties
		NDiff
		NDiffFMO
	end

	methods
		function this = Detector(sz)
			this.Size = sz;
			this.Size_r = round(sz * this.scale_factor);
			this.NDiff = zeros(this.Size_r);
			this.NDiffFMO = zeros(this.Size_r);
		end

		function [] = update_diff(this, delta, fmo)
			this.NDiff = this.NDiff + delta;
			for kk = 1:numel(fmo)
				this.NDiffFMO(fmo(kk).PixelIdxList) = this.NDiffFMO(fmo(kk).PixelIdxList)+1;
			end
		end

		function prob = get_prob(this)
			prob = this.NDiffFMO ./ this.NDiff;
			prob(this.NDiff == 0) = 0.5;
		end

		function loss = diff_loss(im_r,im0_r,regions)
			G_im = imgradient(rgb2gray(im_r));
			G_im0 = imgradient(rgb2gray(im0_r));
			move_im = G_im - G_im0 > 0.2;
			move_im0 = G_im0 - G_im > 0.2;
			% for kkk = 1:2
			% 	move_im = imdilate(move_im, [1 1 1; 1 0 1; 1 1 1]);
			% 	move_im0 = imdilate(move_im0, [1 1 1; 1 0 1; 1 1 1]);
			% end
			loss = zeros(1,numel(regions));
			for k = 1:numel(regions)
				p = sum(move_im(regions(k).PixelIdxList));
				p0 = sum(move_im0(regions(k).PixelIdxList));
				loss(k) = min(p,p0)/regions(k).Length;
			end
		end

		function [fmo, regions_fmo] = detect(this,im, im0, im1)
			sz = [size(im,1) size(im,2)];
			[delta_bin, delta_plus_minus_bin] = this.get_deltas(im,im0,im1);
			regions_fmo = this.get_regions_fmo(delta_bin, delta_plus_minus_bin);
			regions = regions_fmo;

			fmo = [];
			if isempty(regions)
				return;
			end
			dist = bwdist(~delta_bin);
			loss = NaN*ones(size(regions));
			for i = 1:numel(regions)
				regions(i).PixelIdxList = regions(i).PixelIdxList';
				regions(i).PixelList = regions(i).PixelList';
				regions(i).Distances = dist(regions(i).PixelIdxList);
				regions(i).Radius = max(regions(i).Distances);
				normd = regions(i).Distances/regions(i).Radius;
				regions(i).Trajectory = regions(i).PixelIdxList(normd > this.traj_t);
				regions(i).TrajectoryXY = regions(i).PixelList(:,normd > this.traj_t);
				regions(i).Length = 0;
				if ~Frame.is_connected(regions(i).Trajectory,sz)
				   loss(i) = Inf;
				   regions(i).Trajectory = [];
				   regions(i).PixelIdxList = [];
				end
			end
			delta_fmo_bin = logical(zeros(sz));
			delta_fmo_bin([regions.PixelIdxList]) = 1;

			fast = ~isinf(loss);

			bin_traj = logical(zeros(sz));
			bin_traj([regions.Trajectory]) = 1;   
			bin_traj_thin = bwmorph(bin_traj,'thin',Inf);
			for i = 1:numel(regions)
				if isinf(loss(i))
				   continue;
				end
				loss(i) = Inf;

				bin_traj1 = logical(zeros(sz));
				bin_traj1(regions(i).Trajectory) = 1;   
				dist = bwdist(bin_traj1);
				inl = dist(regions(i).PixelIdxList) <= regions(i).Radius;
				
				regions(i).PixelIdxList = regions(i).PixelIdxList(inl);
				regions(i).PixelList = regions(i).PixelList(:,inl);
				regions(i).Distances = regions(i).Distances(inl);

				regions(i).Trajectory = regions(i).Trajectory(bin_traj_thin(regions(i).Trajectory));
				regions(i).Length = numel(regions(i).Trajectory);
			    if regions(i).Length > this.min_length && ...
			       regions(i).Area > this.min_area && ...
			       regions(i).Radius > this.min_radius

				   exp_area = 2*regions(i).Radius*regions(i).Length + pi*regions(i).Radius^2;
				   loss(i) = abs(regions(i).Area/exp_area - 1);
				end
			end
			regions_fmo = regions(fast);
			fmo = regions(loss < this.max_area_dif);
      	end
	end
end
