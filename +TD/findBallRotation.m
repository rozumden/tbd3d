function [u v] = findBallRotation(imgs, timestamps, max_span)
% Image based full search estimation of rotation axis and angular velocity (assumed const for all imags). No interpolations, comparison only in small inner-part of the ball, though the total number of evaluated transforms is large so the processing takes time.
% 
% imgs - cell array of sharp imgs of the ball at different poses during single rotation motion. Can have different sizes but the center of the ball is assumed at the center and the radius is ~min(sz)/2. Must be in the order corredponding to timestamps.
% timestamps - times corresponding to imgs; the return angular velocity will be wrt these units
% max_span - max time difference between poses to be evaluated (rotation should not exceed max_rot); in the same units as timestapms
% 
% u - unit vector of rotation axis orientation; coordinates: x1=positive down, x2=positive right, x3=positive towards you
% v - scalar and always non-negative angular velocity in angle per unit in timestamps (magnitute of the angular velocity vector)

% DBG
% rng(0);

% default args
if(nargin < 2 || isempty(timestamps)) timestamps = 0:length(imgs)-1; end
if(nargin < 3 || isempty(max_span)) max_span = Inf; end
min_eval_r = .5; % percentage of the radius of the "smallest ball" on which the evaluation is calculated - implies upper bound to max rotation between imgs{i} and imgs{j} (see also max_rot). min_eval_r=1 means the whole ball must be visible between poses (-> max_rot=0), while min_eval_r->0 means that only small common part can be visible (estimation can be very inaccurate)
min_span = 0;
max_rot = Inf;

% determine target evaluation sz on which image similarity after rotation will be measured
sz = zeros(length(imgs), 3); % input image 2-size and assumed ball radius
for i=1:length(imgs) sz(i,:) = [size(imgs{i},1) size(imgs{i},2) size(imgs{i},3)]; end
r = (min(sz(:,1:2),[],2)-1)/2; % assumed input image ball radius
offset = (sz(:,1:2)+1)/2; % offset of center-based coords and 1-based matlab coords (basically center position in 1-based coords)
target_sz = min(min(sz(:,1:2),[],2)); target_r = (target_sz-1)/2; % square ball image as output of evaluated transforms

% determine max rotation between imgs
max_rot = min(max_rot, acos(min_eval_r)); % acos - max angle such that the center of the img is visible in both poses
eval_r = floor(cos(max_rot)*target_r);

% setup all evaluated transform combinations (images to transform from->to)
timediff = timestamps(:).'-timestamps(:);
keep = abs(timediff) <= max_span & abs(timediff) > max(min_span,0);
[idx_from idx_to] = find(keep);
from_to = [idx_from idx_to timediff(keep)]; % table of all evaluated transforms: idx of src-dest image and time-difference (positive when past->future)

% discretization of the searched space
v_delta = 2/180*pi; % prescribed samplings spacing for angular velocity search
N_axes = 400; % prescribed number of (full sphere) vectors for axis search (actual tested number is approx N/2, half sphere)

% prepare transformation to normalized pose - rescale to target_sz and extract only indices with the inner circle where transformation error is calculated
temp = -target_r:target_r;
eval_mask = temp(:).^2+temp.^2 <= eval_r^2; % pixels where transformation error will be calculated
[out1 out2] = find(eval_mask); % indices of outer px
out = [temp(out1); temp(out2)].';
out(:,3) = sqrt(target_r.^2-sum(out.^2,2));

% scale all imgs to normalized poses
imgs_normalized = cell(size(imgs));
for i=unique(idx_to).'
	in = round(out(:,1:2)/target_r*r(i) + offset(i,:)); % 1-based indices to imgs{i}
	idx = (in(:,2)-1)*sz(i,1)+in(:,1); idx = idx(:)+(0:sz(i,3)-1)*sz(i,1)*sz(i,2); % linear idx to imgs{i}
	imgs_normalized{i} = imgs{i}(idx); % pixels in evaluated region (reshaped)
end
% interp2 = cat(3, [0 0], [0 1], [1 0], [1 1]);

% best results for each transform
best_u = zeros(size(from_to,1),3); % best unit axis for each transform (scaled such that v>0 when rotation past->future)
best_v = zeros(size(from_to,1),1); % best angular speed for each transform (scaled such that v>0 when rotation past->future)
best_err = Inf(size(from_to,1),1); % actual transformation error of the best combination

% find optimal rotation for each required transforms
for i=1:size(from_to,1)
	ix_from = from_to(i,1); ix_to = from_to(i,2);

	% randomly setup searched space (with similar degree of discretization for all tansforms)
	u = randu(N_axes); % note: all 'u's point towards camera, which is good for most FMO balls captured from the side:-D
	v = mod(linspace(0,2*max_rot,ceil(2*max_rot/v_delta)+1)+rand(1)*2*max_rot, 2*max_rot) - max_rot; % slight randomization of angular velocity sampling

	% try all rotations
	for ix_u = 1:size(u,1)
		for ix_v = 1:length(v)
			u_ = u(ix_u,:); v_ = v(ix_v);
			R = cos(v_)*eye(3)-sin(v_)*[0 -u_(3) u_(2); u_(3) 0 -u_(1); -u_(2) u_(1) 0]+(1-cos(v_))*(u_(:)*u_); % rotation matrix by 'v' (inverse if applied as transpose)

			% transform output->input coords, look up indices
			in = round((out*R(:,1:2))/target_r*r(ix_from) + offset(ix_from,:)); % 1-based indices to imgs{ix_from}; rotation by T.', output z-coord ignored
			% ** interpolated version:
			% in0 = floor(in) + interp2; % 1-based neighboring indices to imgs{ix_from};
			% w = prod(1-abs(in-in0),2); % interpolation weights 
			% in0 = min(max(in0,1),sz(ix_from,[1 2])); % fixes out of range for large rotations; needs to be done after 'w' calculation
			% idx = (in0(:,2,:)-1)*sz(ix_from,1)+in0(:,1,:); idx = idx+(0:sz(ix_from,3)-1)*sz(ix_from,1)*sz(ix_from,2); % linear idx to imgs{ix_from}, incl all RGB channels (dim=2) and interpolation (dim=3)
			% vals = sum(imgs{ix_from}(idx).*w,3); % src image pixels after geometric transform with 2linear interpolation
			% ** non-interpolated version:
			idx = (in(:,2)-1)*sz(ix_from,1)+in(:,1); idx = idx+(0:sz(ix_from,3)-1)*sz(ix_from,1)*sz(ix_from,2); % linear idx to imgs{ix_from}, incl all RGB channels (dim=2) and interpolation (dim=3)
			vals = imgs{ix_from}(idx); % src image pixels after geometric transform

			% scale src pixels to best match dest pixels (y=a*x linear scale, without bias, a\in[.5,1.5])
			scale = (vals(:).'*imgs_normalized{ix_to}(:))/(vals(:).'*vals(:)); % L2 fit
			scale = min(max(scale,.5),1.5); % restrict to sensible range
			%scale = 1; % DBG
			
			% calc error, compare with SOTA
			err = sum(abs(scale*vals(:) - imgs_normalized{ix_to}(:))); % L1 transformation error
			if(err < best_err(i))
				best_u(i,:) = u_; best_v(i) = v_; best_err(i) = err;
			end
		end
	end
end

% normalize 'v' to actual angular velocity
best_v = best_v./from_to(:,3);

% flip axis if v < 0 to align axis vectors
w = best_v < 0;
best_u(w,:) = -best_u(w,:);
best_v(w,:) = -best_v(w,:);

% DBG
% pts = best_u.*best_v/v_gt;
% [X Y Z] = sphere();
% mesh(X,Y,Z, 'facecolor', 'none'); axis equal;
% hold on;
% plot3(pts(:,1), pts(:,2), pts(:,3), 'or');
% plot3(u_gt(:,1), u_gt(:,2), u_gt(:,3), 'oc', 'linewidth', 2);
% hold off;

% transform errors to represent 'score' (ie more is better)
best_err = 1./(best_err+(median(best_err)-min(best_err))/2); % suppresion of poor results without favoring good results too much

score = 0;
inliers = [];
minv = v_delta/(max(timestamps)-min(timestamps)); % minimum discernible rotation? (should be *(num_imgs-1) but this works better...)

% maximum-consensus selection of axis and velocity
for i=1:size(from_to,1)
	u = best_u(i,:); v = best_v(i); % candidate for best estimate, calc consensual score

	% distance of others from the candidate
	d1 = best_u*u(:); % cos dist (angle only)
	d2 = abs(best_v - v)/max(v,minv); % relative err of veliocity only
	%d3 = sum((best_v.*best_u - u*v).^2,2); % euclid dist^2

	% inliers + score
	in = (d1 > .97 & d2 < .1);% | d3 < minv^2; % change inlier criteria here
	s = sum(best_err(in));
	if(s > score)
		score = s;
		inliers = in;
	end
end

% try v=0 as model (does not work well with d1-based metric around zero)
d3 = sum((best_v.*best_u).^2,2); % euclid dist^2 from 0
in = d3 <= minv^2;
s = sum(best_err(in));
if(s > score) score = s; inliers = in; end

% DBG
% fprintf('inliers %d/%d\n', nnz(inliers), size(from_to,1));

% average inliers, return
uv = sum(best_u(inliers,:).*best_v(inliers).*best_err(inliers),1)/score;
v = sqrt(uv*uv(:));
if(v > 1e-3)
	u = uv/v; % unit axis direction

	% average 'v' separately
	v = sum(best_v(inliers).*best_err(inliers))/score;
else
	u = [0 0 1]; v = 0;
end

% DBG
% res = u*v/v_gt;
% hold on;
% plot3(pts(inliers,1), pts(inliers,2), pts(inliers,3), 'xg');
% plot3(res(:,1), res(:,2), res(:,3), 'om', 'linewidth', 2); title(sprintf('cos=%.3f, relv=%.3f', u*u_gt(:), abs(v-v_gt)/v_gt));
% hold off;
end

function x = randu(N)
% evenly distributed and slightly 'randomized' sampling of half-sphere; based on http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
% note:can be 'more random' by additional random rotation
% N: roughly 2x number of points returned (not exact)
% x: cartesian coords of the pts on half-sphere (with x3>0) (roughly N/2 pts, can be N/2+1)

% rnd pts in [0,1]^2
t = (0:N).'./[N (1+sqrt(5))/2];
t = mod(t+rand(1,2),1); % random displacement away from poles
t = t(t(:,1) <= 1/2+.05, :); % keep only half of the sphere (to be) (note: plus little bit more to cover the 'equator')

% map to sphere
ct = 2*sqrt(t(:,1).*(1-t(:,1)));
x = [ct.*cos(2*pi*t(:,2)) ct.*sin(2*pi*t(:,2)) 1-2*t(:,1)];
end
