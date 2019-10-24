function [u rot] = findBallRotation(imgs, timestamps, max_span)
% Image based full search estimation of rotation axis and angular velocity (assumed const for all imags). No interpolations, comparison only in small inner-part of the ball, though the total number of evaluated transforms is large so the processing takes time.
% 
% imgs - cell array of sharp imgs of the ball at different poses during single rotation motion. Can have different sizes but the center of the ball is assumed at the center and the radius is ~min(sz)/2. Must be in the order corredponding to timestamps.
% timestamps - times corresponding to imgs; the return angular velocity will be wrt these units
% max_span - max time difference between poses to be evaluated (rotation should not exceed max_rot); in the same units as timestapms
% 
% u - unit vector of rotation axis orientation; coordinates: x1=positive down, x2=positive right, x3=positive towards you
% rot - scalar and always positive angular velocity in angle per unit in timestamps (magnitute of the angular velocity vector)

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

% setup combinations of images to tronsform from-to
timediff = timestamps(:).'-timestamps(:);
keep = abs(timediff) <= max_span & abs(timediff) > max(min_span,0);
[idx_from idx_to] = find(keep);
from_to = [idx_from idx_to timediff(keep)]; % table of all evaluated transforms: idx of src-dest image and time-difference (positive when past->future)
idx_from = unique(idx_from(:)).'; idx_to = unique(idx_to).';

% discretization of the searched space
% rot_delta = atan(.5/target_r); % difference of 1/2px at max radius
rot_delta = 2/180*pi;
rot = linspace(0,max_rot,ceil(max_rot/rot_delta)); rot = [-rot(end:-1:2) rot];
tilt_dir = (0:10:350)/180*pi; % equidistant spacing for axis orientation
tilt_amount = (0:10:90)/180*pi; % equidistant spacing for axis inclination
err = NaN(length(tilt_dir), length(tilt_amount), length(rot), size(from_to,1)); % transformation error from i-th to j-th img. Indexed as (tilt_dir, tilt_amount, rot angle, transform_idx)
do_search = false(size2(err,1:3)); % mask of axis orientation and rotatons which will be evaluated
do_search(1,1,:) = 1; % dir=0 for tilt = 0
do_search(1:6:end,2,:) = 1; % dir=0 60 120 ... for tilt = 10
do_search(1:3:end,3,:) = 1; % dir=0 30 60 ... for tilt = 20
do_search(1:2:end,4,:) = 1; % dir=0 20 40 ... for tilt = 30
do_search(:,5:end-1,:) = 1; % all dirs for tilt >= 40 < 90
do_search(1:18, end, :) = 1; % half circle for tilt=90 (axis can be flipped)
do_search(:,:,rot == 0) = 0; do_search(1,1,rot == 0) = 1; % only 1 axis orientation for no rotation

% prepare transformation to normalized pose - rescale to target_sz and extract only indices with the inner circle where transformation error is calculated
temp = -target_r:target_r;
eval_mask = temp(:).^2+temp.^2 <= eval_r^2; % pixels where transformation error will be calculated
[out1 out2] = find(eval_mask); % indices of outer px
out = [temp(out1); temp(out2)].';
out(:,3) = sqrt(target_r.^2-sum(out.^2,2));

% scale all imgs to normalized poses
imgs_normalized = cell(size(imgs));
for i=idx_to
	in = round(out(:,1:2)/target_r*r(i) + offset(i,:)); % 1-based indices to imgs{i}
	idx = (in(:,2)-1)*sz(i,1)+in(:,1); idx = idx(:)+(0:sz(i,3)-1)*sz(i,1)*sz(i,2); % linear idx to imgs{i}
	imgs_normalized{i} = imgs{i}(idx); % pixels in evaluated region (reshaped)
end

% search all axis orientations
for ix=find(do_search).'
	[ix_dir ix_amount ix_rot] = ind2sub(size2(err,1:3), ix);

	% rotation matrix (coordinates: x1=positive down, x2=positive right, x3=positive towards you)
	t1 = tilt_dir(ix_dir); t2 = tilt_amount(ix_amount); % axis tilt
	u = [-cos(t2) sin(t2)*sin(t1) sin(t2)*cos(t1)]; % unit vector of the axis of rotation (in x1 x2 x3 coords)
	T = cos(rot(ix_rot))*eye(3)-sin(rot(ix_rot))*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(rot(ix_rot)))*(u(:)*u); % rotation matrix by 'rot' (inverse if applied as transpose)

	% transform output->input coords (without scaling)
	in_small = out*T(:,1:2); % rotation by T.', output z-coord ignored
	% assert(all(in_small(:,3) >= 0), 'chybka - negative z po rotaci'); % DBG

	% transform each img to every other, calc error
	for i=idx_from
		% scale coords, look up indices
		in = round(in_small/target_r*r(i) + offset(i,:)); % 1-based indices to imgs{i}
		idx = (in(:,2)-1)*sz(i,1)+in(:,1); idx = idx(:)+(0:sz(i,3)-1)*sz(i,1)*sz(i,2); % linear idx to imgs{i}
		vals = imgs{i}(idx); % i-th image pixels after transform

		% compare agains every other img
		for j=find(from_to(:,1) == i).'

			% scale src pixels to best match dest pixels (y=a*x linear scale, without bias, a\in[.5,1.5])
			scale = (vals(:).'*imgs_normalized{from_to(j,2)}(:))/(vals(:).'*vals(:)); % L2 fit
			scale = min(max(scale,.5),1.5); % restrict to sensible range
			
			% calc error
			err(ix_dir, ix_amount, ix_rot, j) = sum(abs(scale*vals(:) - imgs_normalized{from_to(j,2)}(:))); % L1 transformation error
		end
	end
end

% aggregate common axis orientation for all transforms and find best
err_tilt = sum(min(err, [], 3), 4); % for each axis choose best rotation angle of each transform, sum over all img transforms (results in 2d array)
temp = find(any(do_search,3));
[~, t] = min(err_tilt(temp)); [t1 t2] = ind2sub(size(err_tilt), temp(t));

% for each transformation pair choose the best rotation angle for the given axis
[~, r] = min(err(t1, t2, :, :), [], 3); r = squeeze(r);

% 'robustly' estimate single angular velocity for all imgs
x = from_to(:,3).'; y = rot(r); % will fit y=a*x where a is the angular velocity
keep = abs(y) < max_rot-1e-3 & abs(y) > 1e-3; % discard too big and too small values
% x0 = x; y0 = y; % DBG
x = x(keep); y = y(keep);

% initial ransac line fit to avoid outliers
if(~isempty(x))
	keep = ransac_scale(x,y,1.5*rot_delta);
	x = x(keep); y = y(keep);

	% (roughly) L1 line fit to x->y to find linear coeff of the angular velocity
	if(~isscalar(x))
		rot = L1_scale(x,y);
	else
		if(x == 0) rot = 0; else rot = y./x; end
	end
	
	% DBG
	% figure(1); plot(x0,y0, x0, rot*x0, x, y, 'og', x0, rot*x0+2*rot_delta, '--k', x0, rot*x0-2*rot_delta, '--k'); title(sprintf('dir=%d, tilt=%d, rot=%.2f deg/time1', round(tilt_dir(t1)/pi*180), round(tilt_amount(t2)/pi*180), rot/pi*180));
	% drawnow;
else
	rot = 0;
end

% recalc orientation to unit axis angle
u = [-cos(tilt_amount(t2)) sin(tilt_amount(t2))*sin(tilt_dir(t1)) sin(tilt_amount(t2))*cos(tilt_dir(t1))]; % unit vector of the axis of rotation (in x1 x2 x3 coords); coordinates: x1=positive down, x2=positive right, x3=positive towards you

% flip axis if rot < 0 for easier aligning and filtering
if(rot < 0)
	rot = -rot;
	u = -u;
end
end

function inliers = ransac_scale(x,y,inlier_thresh)
% fits y=a*x model and identifies inliers

% remove x==0 - always inlier and causes problems with /0
inliers = false(size(x));

% MSS is 1 point, do full search
for i=1:length(x)
	if(x(i) == 0) continue; end
	in = abs((y(i)*x)/x(i)-y) <= inlier_thresh; % approx error and inliers
	if(nnz(in) > nnz(inliers))
		inliers = in;
	end
end
end

function a = L1_scale(x,y)
% roughly L1 fit y=a*x

w = 1; % L1 reweighting
a = 0; % scaling factor
for iter = 1:5 % max iter
	a_old = a;
	wx = w.*x;
	a = wx*y(:)/(wx*x(:));
	if(abs(a_old-a)/a < 1e-3)
		break;
	end
	w = 1./sqrt((a.*x(:)-y(:)).^2+1e-8).';
end
end
