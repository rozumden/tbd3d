function [coeffs img val fits models psf_mask log] = psffit(psf,params)
is_dbg = false;
log = string([]);
thresh1 = 55/255; thresh2 = 30/255; % linear
gamma = 2;
quad_sampl_thresh = 40/255; % value (linear, relative to segment max)
quad_max_pts = 40; % max number of pts to sample from for quadratic ransac (chooses approx this many best pts at most)
max_quad_curvature = .200;
min_quad_length = 10;
val_disk_sz = 5; % diameter

% lin ransac params
sig_thresh1 = (80/255)^gamma; sig_thresh2 = (45/255)^gamma; % linear->gamma
gap_thresh_small = 1.5; gap_thresh_big = 5+1e-6;
sig_count1 = 4; sig_pcent1 = .4; % ^gamma
sig_count2 = 6; sig_pcent2 = .8; % ^gamma

% % quad ransac params
% sig_thresh1 = (80/255)^gamma; sig_thresh2 = (45/255)^gamma; % linear->gamma
% gap_thresh_small = 1.5; gap_thresh_big = 5+1e-6;
% sig_count1 = 4; sig_pcent1 = .4; % ^gamma
% sig_count2 = 6; sig_pcent2 = .8; % ^gamma

% fit candidate collection
fits = {};

% threshold psf
psf(psf < 0) = 0; psf = psf/max(psf(:)); % normalize to [0,1]
psf_mask = psf_supp(psf >= thresh1, ones(3,3), psf >= thresh2); % psf 'support' which will be fitted
[x1 x2] = find(psf_mask); pts = [x1 x2]; % fitted pts
psf_gamma = psf.^gamma;
w = psf(psf_mask); w_gamma = psf_gamma(psf_mask);
se = diskMask(val_disk_sz); % strel for determining fit value (coverage of PSF)

% find all lin segments by ransac
if(is_dbg) rng(0); end % DBG - deterministic behavior
models = seqransac_lin(pts, w_gamma, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, params);

if(~isempty(models))
	% re-fit best lin model
	control_pts = models(1).end_pts;
	[coeffs metrics] = fitpwpoly_ng(pts, control_pts, 1, w_gamma, 200, 1, 2);

	% trim (remove small insignificant garbage on the sides)
	if(~isempty(coeffs))
		inliers = find(metrics.inliers);
		[cs] = findSignificantRun(metrics.curve_s(inliers), w_gamma(inliers), ones(numel(inliers),1), [], gap_thresh_small, 0, sig_thresh1, sig_count1, 0, sig_thresh2, sig_count2, 0);
		if(~isempty(cs))
			t0 = metrics.curve_t(inliers(cs(1))); t1 = metrics.curve_t(inliers(cs(end))); % parametrization of the new curve-ends (cropped)
			if(t1 > t0)
				c = coeffs{1}; % old coeffs
				coeffs = {[c(1,:)+c(2,:)*t0; c(2,:)*(t1-t0)]}; % reparametrization [t0,t1]->[0,1]

				% add to fit-candidate collection
				fits{end+1} = struct('id', 'lin1', 'coeffs', {coeffs}, 'control_pts', control_pts);
			elseif(is_dbg)
				log(end+1) = 'lin1 trim->dot';
			end
		elseif(is_dbg)
			log(end+1) = 'lin1 trim->empty';
		end
	elseif(is_dbg)
		log(end+1) = 'lin1 fitpwpoly failed';
	end

	% quadratic
	if(~isscalar(models)) % use more segments as base for quad sampling
		idx1 = models(1).cs; len1 = models(1).len; % primary segment cs and len
		idx2 = [models(2:min(3,end)).cs]; % take at moset 3 best models to consideration
		len2 = max([models(2:min(3,end)).len]); % length of the 'secondary' segments ("segment2")

		
		count1 = len1*quad_max_pts/(len1+len2); count2 = ceil(quad_max_pts-count1); % number of sampling pts used from model1/2 such in the ratio of the respective lengths
		temp = w(idx2); idx2 = idx2(temp > quad_sampl_thresh*max(temp(:))); % sufficient intensity
		idx2 = idx2(quad_find_ransac_buddies(models(1), pts(idx2,:), max_quad_curvature/2, 2)); % suitable position wrt first segment; curv/2 - be even more strict in the initial sampling
		if(numel(idx2) > count2)
			if(is_dbg) rng(0); end % DBG - deterministic behavior
			idx2 = idx2(randperm(numel(idx2), count2)); % randomly select small subset
		end
		
		% get pts from the primary segment1
		count1 = quad_max_pts-numel(idx2);
		temp = w(idx1); idx1 = idx1(temp > quad_sampl_thresh*max(temp(:))); % sufficient intensity
		if(numel(idx1) > count1)
			if(is_dbg) rng(0); end % DBG - deterministic behavior
			idx1 = idx1(randperm(numel(idx1), count1)); % randomly select small subset
		end
		idx = [idx1 idx2];
	else % use 1st segment only
		idx = models(1).cs;
		temp = w(idx); idx = idx(temp > quad_sampl_thresh*max(temp(:))); % sufficient intensity
		if(numel(idx) > quad_max_pts)
			if(is_dbg) rng(0); end % DBG - deterministic behavior
			idx = idx(randperm(numel(idx), quad_max_pts)); % randomly select small subset
		end
	end
	if(numel(idx) >= 4)
		if(is_dbg) rng(0); end % DBG - deterministic behavior
		model_quad = ransac_quad(pts, w_gamma, idx, max_quad_curvature, gap_thresh_small*1.3, gap_thresh_big*1.5, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, params); % dist thresholds slightly increased wrt to linear

% 		% DBG
% 		temp = find(psf_mask); m = false(size(psf_mask)); m(temp(idx)) = true;
% 		imshow(imoverlay(psf, m, [1 0 0], .4));
% 		if(model_quad.val)
% 			hold on;
% 			p = linspace(model_quad.t_bounds(1), model_quad.t_bounds(2), 20).'.^(0:2)*model_quad.coeffs;
% 			plot(p(:,2), p(:,1), '-');
% 			cp = ([1 0; .5 .5; 0 1]*model_quad.t_bounds).^(0:2)*model_quad.coeffs;
% 			plot(cp(:,2), cp(:,1), 'o');
% 			plot(model_quad.mss(:,2), model_quad.mss(:,1), 'g+');
% 			hold off;
% 		end
% 		pause;

		if(model_quad.val > 0 && model_quad.len >= min_quad_length)
			control_pts = ([1 0; .5 .5; 0 1]*model_quad.t_bounds).^(0:2)*model_quad.coeffs; % start-mid-end of the parabolic arc
			[coeffs metrics] = fitpwpoly_ng(pts, control_pts, 2, w_gamma, 200, 1, 2, [80 6]);

			if(~isempty(coeffs))
				% check final metrics
				if(metrics.length >= min_quad_length && metrics.max_curvature <= max_quad_curvature) % note: length should be measured after trimming but this is simpler...
					% trim (remove small insignificant garbage on the sides)
					inliers = find(metrics.inliers);
					[cs] = findSignificantRun(metrics.curve_s(inliers), w_gamma(inliers), ones(numel(inliers),1), [], gap_thresh_small*1.3, 0, sig_thresh1, sig_count1, 0, sig_thresh2, sig_count2, 0); % dist thresholds slightly increased wrt to linear
					if(~isempty(cs))
						t0 = metrics.curve_t(inliers(cs(1))); t1 = metrics.curve_t(inliers(cs(end))); % parametrization of the new curve-ends (cropped)
						if(t1 > t0)
							c = coeffs{1}; % old coeffs
							coeffs = {[t0.^(0:2)*c; (t1-t0)*(c(2,:)+2*t0*c(3,:)); (t1-t0)^2*c(3,:)]}; % reparametrization [t0,t1]->[0,1]

							% add to fit-candidate collection
							fits{end+1} = struct('id', 'quad1', 'coeffs', {coeffs}, 'control_pts', control_pts, 'max_curvature', max(metrics.max_curvature));
						elseif(is_dbg)
							log(end+1) = 'quad1 trim->dot';
						end
					elseif(is_dbg)
						log(end+1) = 'quad1 trim->empty';
					end
				elseif(is_dbg)
					log(end+1) = sprintf('quad1 post-rejected (len=%.1f, curv=%.2f)', metrics.length, metrics.max_curvature);
				end
			elseif(is_dbg)
				log(end+1) = 'quad1 fitpwpoly failed';
			end
		elseif(is_dbg)
			if(model_quad.val)
				log(end+1) = sprintf('quad1 pre-rejected (len=%.1f)', model_quad.len);
			else
				log(end+1) = sprintf('quad1 not found');
			end
		end
	elseif(is_dbg)
		log(end+1) = sprintf('quad1 skipped (#pts=%d)', numel(idx));
	end
end
if(isempty(fits))
	% no good fit possible, return delta
	sz = 3; % velikost okoli
	cc = conv2(psf_gamma, double(diskMask(sz)), 'valid');
	[x1 x2] = find(cc == max(cc(:))); % best delta placement absed on val_gamma (topleft)
	pos = [x1 x2] + (sz-1)/2; % delta pos (possibly several maxima)
	[~,ix] = max(psf((pos(:,2)-1)*size(psf,1)+pos(:,1))); % find overall best placement
	fits{end+1} = struct('id', 'delta', 'coeffs', {{pos(ix,:)}}, 'control_pts', []);
end

% determine goodness of fit
for i=1:length(fits)
	[fits{i}.img fits{i}.val fits{i}.valg fits{i}.bleed] = fit_val(fits{i}.coeffs, se, psf, psf_gamma, psf_mask);
end

% heuristic decision between lin and quad
if(length(fits) > 1)
	%if((fits{1}.val+.05 > fits{2}.val && fits{2}.max_curvature < .015) || (fits{1}.val+.10 > fits{2}.val && fits{2}.max_curvature > .100))
	if((fits{1}.val+.10 > fits{2}.val && fits{2}.max_curvature > .100))
		ix = 1; % prefer lin
	else
		[~,ix] = max([fits{1}.val, fits{2}.val]);
	end
else
	ix = 1;
end
coeffs = fits{ix}.coeffs;
img = fits{ix}.img;
val = fits{ix}.val;
end

function supp = psf_supp(base, dilate, other)
labels = bwlabel(other);
supp = false(size(base));
temp = unique(labels(imdilate(base, dilate))).';
for i=temp(temp > 0)
	supp = supp | labels == i;
end
end

function keep = quad_find_ransac_buddies(m1, pts, max_curvature, bandwidth)
% experimental - finds which points among pts are 'suitable' as ransac-seeds for quadratic fit based on m1 to get parabola of bounded curvature (only very approximate and too conservative estimate - more points can be connected by a parabola than what is returned)
% The allowed region looks a bit like magnet fields lines (protrudng from the ends and not much to the side), but parabolic (not closed):-D
% This criterion fails eg in case of intersecting segments, parts of which together form a parabola, or laterally far-away segments oriented in such a way that a parabola is possible (both of these cases will be rejected)
% 
% bandwidth - sth like inlier_thersh - half width of the band alond m1 and offset of the the parabolas at the ends

% project end_pts onto the first segment
t = (pts-m1.end_pts(1,:))*m1.d(:); % t-param val of the projection (actual length in px)
n = abs((pts-m1.end_pts(1,:))*[m1.d(2); -m1.d(1)]); % perpendicular distance from line of m1

% three cases - pt is 'left' of m1 (t<0), pt is 'above' m1 (t in [0,len]), and pt is 'right' of m1 (t>len)
set1 = t < 0; set3 = t > m1.len; set2 = ~(set1 | set3);
keep = false(size(pts,1),1);
keep(set1) = max_curvature*t(set1).^2+bandwidth >= n(set1);
keep(set3) = max_curvature*(t(set3)-m1.len).^2+bandwidth >= n(set3);
keep(set2) = bandwidth >= n(set2);
end

function [img val valg bleed] = fit_val(coeffs, se, psf, psf_gamma, psf_mask)
img = logical(renderpwpoly2(coeffs, size(psf)));
cover = imdilate(img, se) & psf_mask;
val = sum(psf(cover))/sum(psf(psf_mask));
valg = sum(psf_gamma(cover))/sum(psf_gamma(psf_mask));
bleed = nnz(img & ~psf_mask);
end
