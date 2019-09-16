% Abstract class for creating FMO detector
classdef (Abstract) IDetector < handle
   properties
      noise_t_b = 10/255
      min_noise_t_b = 4/255
      max_noise_t_b = 12/255
      noise_inc = 2/255

      noise_t = [0.1 0.08 0.11]/2
      noise_t_addapted = []
      fmo_t = 0.9
      traj_t = 0.7

      min_area = 10
      noise_area = 3
      min_radius = 2
      min_length = 10
      max_area_dif = 0.6
      strokiness_th = 0.7

      do_stabilize = false
      do_post_process = true
      scale_factor = []
      max_height = 960

      Size
      Size_r
   end

   methods (Abstract)
      detect(this, im)
   end
   
   methods
      function bin_out = binarize_addaptive(this, delta)
         acc = sum(delta,3);
         while true
            bin = acc > this.noise_t_b;
            bin_no_noise = bwareaopen(bin,this.noise_area);
            bin_out = bin_no_noise;
            noise = bin & ~bin_no_noise;
            if this.noise_t_b > this.min_noise_t_b & sum(noise(:)) < 70 
               this.noise_t_b = this.noise_t_b - this.noise_inc;
            elseif this.noise_t_b < this.max_noise_t_b & sum(noise(:)) > 150 
               this.noise_t_b = this.noise_t_b + this.noise_inc;
               break;
            else
               break;
            end
         end
         sortedNoise = sort(acc(noise));
         if isempty(sortedNoise)
            this.noise_t_addapted = this.noise_t_b;
         else
            this.noise_t_addapted = sortedNoise(ceil(0.8*end));
            for k = 1:2
               dt = bwdist(~bin_out);
               bin_out(dt == 1 & acc < this.noise_t_addapted) = 0;
            end
            bin_out = bwareaopen(bin_out,this.noise_area);
         end
      end

      function bin = binarize_stochastic(this, differ)
         acc = sum(differ,3);
         [a,b] = imhist(acc); % P(differ = theta)
         probs = cumsum(a)/sum(a);
         ind = find(probs > 0.9);
         theta_small = b(ind(1));

         [~,m_ind] = max(a);
         theta_max = 4*b(m_ind);

         speed_of_change = (a - circshift(a,-2))./(a+eps);
         speed_of_change(1:m_ind) = 1; 
         ind = find(abs(speed_of_change) < 0.15);
         theta_flat = b(ind(1));
         theta = max([1/255 theta_max theta_flat theta_small]);
         bin = acc > theta;
         for k = 1:1
            dt = bwdist(~bin);
            bin(dt < 1.5 & acc < 2*theta) = 0;
         end
         for k = 1:1
            dt = bwdist(bin);
            bin(dt < 1.5 & acc > theta/2) = 1;
         end
         bin = bwareaopen(bin,this.noise_area);
         % bin = ~bwareaopen(~bin,this.noise_area);
         this.noise_t_b = theta;
      end


      function bin = binarize(this, delta)
         bin = sum(delta,3) > this.noise_t_b;
         % bin = delta(:,:,1) > this.noise_t(1) | ...
         %       delta(:,:,2) > this.noise_t(2) | ...
         %       delta(:,:,3) > this.noise_t(3);
         % bin = rgb2gray(delta) > this.noise_t_b;
      end

      function bin = binarize_strict(this, delta)
         acc = sum(delta,3);
         bin = acc > 3*this.noise_t_b;
         for k = 1:1
            dt = bwdist(~bin);
            bin(dt == 1 & acc < 6*this.noise_t_b) = 0;
         end
      end

      function [delta_bin, delta_pm_bin, delta_plus_bin, delta_minus_bin] = get_deltas(this,im,im0,im1)
         if this.do_stabilize
            im1t = stabilize(im,im1);
            im0t = stabilize(im,im0);
         else
            im1t = im1;
            im0t = im0;
         end

         delta_plus = abs(im - im0t);
         delta_minus = abs(im1t - im);
         delta0 = abs(im1t - im0t);

         delta_plus_bin = this.binarize(delta_plus);
         delta_minus_bin = this.binarize(delta_minus);
         delta0_bin = this.binarize(delta0);

         delta_pm_bin = delta_plus_bin & delta_minus_bin;
         if this.do_post_process
            delta_pm_bin = Differential.post_process(delta_pm_bin);
         end

         delta_bin = delta_pm_bin & ~delta0_bin;
         if this.do_post_process
            delta_bin = Differential.post_process(delta_bin);
         end
      end

      function regions_fmo = get_regions_fmo(this,delta_bin, delta_pm_bin)
         regions_fmo = regionprops(delta_pm_bin, 'MajorAxisLength', ...
            'PixelList', 'PixelIdxList','Area', 'Orientation', 'BoundingBox');
         scores = arrayfun(@(x) sum(delta_bin(x.PixelIdxList))/x.Area,regions_fmo);
         regions_fmo = regions_fmo(scores > this.fmo_t);
      end

   end
end 