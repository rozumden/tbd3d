function model = ransac_quad(pts, w, idx2sample, max_curvature, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, params)
% idx2sample - indices to 'pts', which consider as mss

% params
eps = 1e-6;
inlier_thresh = 2+eps;
% gap_thresh = 6+eps; % note: gap=1 (or sqrt(2) in euclidian distance) is actually no gap in the pixel grid
max_samples = params.ransac_quad_max_samples;

% generate samples from idx2sample
num_samples = nchoosek(numel(idx2sample), 4); % total number of possible quadruples to sample
if(num_samples > max_samples)
	idx = randsample(num_samples, max_samples);
else
	idx = 1:num_samples; % try all
end
idx = nthtuple(idx, 4); % indices to idx2sample

% component labeling (connectivity between pieces of psf)
sz = max(pts,[],1)+2; % psf size (cropped)
temp = accumarray(pts, 1, sz); temp = temp | imclose(temp, ones(3)); % psf super-mask
[labels num] = bwlabel(temp); labels = labels((pts(:,2)-1)*sz(1)+pts(:,1)); % component label of each pt
lbl_sum = accumarray(labels, w, [num 1]); % value of individual components

% look for max CS
model.val = 0; % currently best model
for i = 1:size(idx,1)
	m = findCS_arc(pts(idx2sample(idx(i,:)), :), pts, w, labels, lbl_sum, inlier_thresh, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, max_curvature, model.val);
	if(m.val > model.val)
		model = m; % currently best model
	end
end
end

function model = findCS_arc(mss, pts, w, labels, lbl_sum, inlier_thresh, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, max_curvature, min_val)

model.val = 0; % bust

% find the two parabolas corresponding to the mss
coeffs_implicit = parabola4p(mss);
if(isempty(coeffs_implicit)) return; end % no real non-degenerate parabola

% heuristic - choose the flatter of the two parabolas (the one where the arc spanned by mss has lower max curvature) - usually one of the parabolas is 'reasonable' and the other 'wrong' one is very narrow
curvature = zeros(size(coeffs_implicit,1),1);
coeffs_param = cell(1,size(coeffs_implicit,1)); R = cell(1,size(coeffs_implicit,1)); coeffs = cell(1,size(coeffs_implicit,1));
for i=1:size(coeffs_implicit,1)
	[coeffs_param{i} R{i} coeffs{i}] = implicit2parametric(coeffs_implicit(i,:)); % coeffs - parabola in the (isometrically) rotated 'canonic' space y~x^2
	% mss_x = mss*R{i}(1,:).'; % x coords of the mss in the rotated space
	% x = min(max(-coeffs{i}(2)/(2*coeffs{i}(3)), min(mss_x)), max(mss_x)); % the pt in arc spanned by mss where max curvature is attained
	x = -coeffs{i}(2)/(2*coeffs{i}(3)); % pt where globally max curvature is attained (different version)
	curvature(i) = 2*abs(coeffs{i}(3))/(1+(2*coeffs{i}(3)*x+coeffs{i}(2))^2)^(3/2); % max curvature
end

[curvature ix] = min(curvature);
if(curvature > 2*max_curvature) % some tolerance allowed, will be screened later
	return;
end
coeffs_param = coeffs_param{ix}; R = R{ix}; coeffs = coeffs{ix}; coeffs_implicit = coeffs_implicit(ix,:);

% % DBG
% plot(pts(:,2), pts(:,1), '+k'); axis equal;
% hold on;
% plot(mss(:,2), mss(:,1), 'ok');
% dbg_plot_implicit_parabola(coeffs_implicit, [max(pts(:,1)) max(pts(:,2))]+10, 'b');
% hold off;

% find all pts that are inliers to the unbounded parabola
[cs t s] = findCS_full(pts*R.', coeffs, inlier_thresh); % rotate pts to the canonical pose where 'coeffs' is regular high-school coefficients of the parabola ~ y=x^2

% find max-value consecutive arc
[cs2 val] = findSignificantRun(s, w(cs), labels(cs), lbl_sum, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2);

if(val > min_val)
	model.val = val;
	model.len = s(cs2(end))-s(cs2(1)); % length along the parabolic arc
	model.cs = cs(cs2); % idx to pt, consensus set
	model.coeffs = coeffs_param; % parametric coeffs (as in fitpwpoly etc)
	model.t_bounds = t(cs2([1; end])); % interval of parameter values correspoinding to the arc alone - the coeffs can be reparametrized to [0,1] later if necessary
	model.coeffs_implicit = coeffs_implicit; % dbg
	model.mss = mss; % dbg
end
end

function [cs t_val s_val] = findCS_full(pts, coeffs, inlier_thresh)
% Finds csonsensus-set for the given parabola - points which lie closer than inlier_thresh to the parabola - and their parametrization and distance along the parabola
%
% cs - indices to pts, consensus set (inliers withing the threshold tolerance)
% t_val - parametrization of pts in cs wrt coeffs_param (not in [0,1] or anything, can be arbitrary values)
% s_val - distance of pts in cs along the parabola (arc-length) from some unspecified zero point (intended for relative mutual distances)

maxiter = 10;
abs_tol = 1e-2; % between-iter stopping condition; abs displacement in either x- or y- of the projected point (max over all pts)

% horizontal projection of each input pt
% NOTE: bounding projections to x_min and x_max in each iteration (which requires precomputing the horizontal projections etc) is not any faster a the initial overhead probably doesn't pay off (as well is the initialization with the horizontal projection). Maybey it is 'safer' that the method doesn't diverge but it probably isn't faster. So if speed is an issued, this should be tested further.
vertex = [-coeffs(2)/(2*coeffs(3))  -coeffs(2)^2/(4*coeffs(3))+coeffs(1)]; % xy coord of the parabola vertex (axis and max)
x2 = zeros(size(pts,1),1); x2(:) = vertex(1); % xy-coordinate of the horizontal projection of each pt (populated later, parabola vertex is the default for pts without horizontal proj)
%y2 = zeros(size(pts,1),1); y2(:) = vertex(2); % xy-coordinate of the horizontal projection of each pt (populated later, parabola vertex is the default for pts without horizontal proj)

discr = coeffs(2)^2-4*coeffs(3)*(coeffs(1)-pts(:,2)); % diskriminant of the eq for the horizontal projection
where = discr > 0;
is_right = pts(:,1) > vertex(1); % pt on the 'right' side of the parabola
factor = (2*is_right(where)-1)*sign(coeffs(3)); % +1 for pts right of the axis, -1 for pts left of the axis; then negated if the parabola is upside down (this says direction of the 'closer' of the two horizontal projections for each pt)
x2(where) = (-coeffs(2)+factor.*sqrt(discr(where)))/(2*coeffs(3));
%y2(where) = pts(where, 2); % actual horizontal projection of each input pt (the nearer of the two solutions, due to 'factor')

% determine bounds of the projection in 'x'
x_min = pts(:,1); x_max = x2;
swap = x_min > x_max;
x_min(swap) = x_max(swap); % bounds of the projection of each pt (the true proj is in between)
x_max(swap) = pts(swap, 1); % bounds of the projection of each pt (the true proj is in between)

% % ### DBG ###
% plot_bounds = [min(pts(:,1)) max(pts(:,1))]; plot_bounds = plot_bounds+[-1 1]*.1*(plot_bounds(2)-plot_bounds(1)); % plot range
% temp = linspace(plot_bounds(1), plot_bounds(2), 100);
% plot(temp, temp(:).^(0:2)*coeffs(:), '-b'); axis equal;
% hold on;
% %y_min = x_min.^(0:2)*coeffs(:); y_max = x_max.^(0:2)*coeffs(:);
% %plot([pts(:,1) x_min].', [pts(:,2) y_min].', '-r');
% %plot([pts(:,1) x_max].', [pts(:,2) y_max].', '-g');
% y = x.^(0:2)*coeffs(:);
% plot([pts(:,1) x].', [pts(:,2) y].', '-r');
% hold off;
% return;

% manhattan distance to parabola
xy = [pts(:,1) pts(:,1).^(0:2)*coeffs(:)]; % actual vertical projection of input (x-coord will be used later, this is also preallocation for the initialization:-)
dist1 = abs(pts(:,2) - xy(:,2)); % vertical distance
%dist2 = abs(pts(:,1) - x2); % horizontal distance

% heuristic - discard some pts based on their algebraic distance
dmax = x_min; dmax(is_right) = x_max(is_right); % x-coord of the maximum derivative for the given pt in its range given by x_min->x_max
dmax = 2*coeffs(3)*dmax + coeffs(2); % actual max derivative
keep = dist1 < inlier_thresh*sqrt(dmax.^2+1); % excluded pts have algebraic distance such that the geometric distance cannot fall below inlier_thresh
if(coeffs(3) > 0)
	keep = keep & pts(:,2) > vertex(2)-inlier_thresh; % these pts lie way bellow parabola (below vertex level)
else
	keep = keep & pts(:,2) < vertex(2)+inlier_thresh; % these pts lie way bellow parabola (below vertex level)
end
%hold on;plot(pts(~keep,2), pts(~keep,1), 'b+');hold off; % DBG
cs = find(keep); % current cs candidates
xy = xy(cs, :); pts = pts(cs, :); x_min = x_min(cs); x_max = x_max(cs); % discard pts and auxiliaries
% fprintf('%d/%d pts rejected upfront\n', numel(keep)-nnz(keep), numel(keep)); % DBG

% switch initialization to horizontal projection if it is closer
%where = dist2(keep) < dist1(keep);
%xy(where,:) = [x2(where) y2(where)];

% iterative projection
for iter = 1:maxiter
	xy_old = xy;
	% % ### DBG - intermediate plot ###
	% temp = linspace(plot_bounds(1), plot_bounds(2), 100);
	% plot(temp, temp(:).^(0:2)*coeffs(:), '-b'); axis equal; axis([min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))] + 60*[-1 1 -1 1]); title(sprintf('after %d iters', iter-1)); %axis([0 100 0 100]);
	% hold on;
	% plot([pts(:,1) xy(:,1)].', [pts(:,2) xy(:,2)].', '-+r');
	% hold off;
	% pause;
	
	% step
	dy = 2*coeffs(3)*xy(:,1)+coeffs(2); % derivative at the current projection estimate
	d = xy - pts; % projection - pt
	xy(:,1) = xy(:,1) - (d(:,1) + dy.*d(:,2))./(dy.^2 + 1 + 2*max(coeffs(3)*d(:,2),0)); % hybrid step (a mix of gradient-descent and newton - takes 'smaller' of the two steps (not literally))
	%xy(:,1) = xy(:,1) - (d(:,1) + dy.*d(:,2))./(dy.^2 + 1 + 2*coeffs(3)*d(:,2)); % pure newton
	
	% clipping
	xy(:,1) = min(max(xy(:,1), x_min), x_max);
	
	% eval
	xy(:,2) = xy(:,1).^(0:2)*coeffs(:); % parabola pt for the current projection estimate

	% exit
	diff = max(abs(xy_old(:)-xy(:)));
	% fprintf('iter %d, diff=%.3f\n', iter, diff); % DBG
	if(diff < abs_tol) break; end
end

% distance thresholding - final CS
keep = sum((pts - xy).^2,2) <= inlier_thresh^2;
cs = cs(keep);

% parameter values for the pts in the final CS
t_val = xy(keep,1);

% distance values along the parabola (basically the inverse of arc-length parametrization)
if(abs(coeffs(3)) > 1e-10) % not a straight line
	temp = 2*coeffs(3)*t_val+coeffs(2);
	temp2 = sqrt(temp.^2+1);
	s_val = (log(temp+temp2)+temp.*temp2)/(4*coeffs(3)); % follows from arc-length integration
else % straight-line approximation
	s_val = sqrt(1+coeffs(2)^2)*t_val;
end
end

function [coeffs_parametric R coeffs_new] = implicit2parametric(coeffs)
% returns parametic repezentation [x y] = sum_i( coeffs_parametric(i)*t^i ) of a parabola given by implicit 'coeffs' (see parabola4p for coeffs format)
% 
% coeffs - coeffs of the implicit parabola f(x,y)=0 as returned by parabola4p
% coeffs_parametric - standard parametric coeffs as used eg by fitpwpoly etc
% R - rotation (orthogonal) matrix transforming the original parabola space to 'canonic' one where y~x^2 (R.' transforms new->old)
% coeffs_new - simple parabola coeffs in the new rotated space of the form y=c2^2+c1*x+c0

% non-scaled version - implicit coeffs do not need to be prescaled
scale2 = sum(coeffs(1:2).^2); scale = sqrt(scale2); % the input coeffs can be scale so that a^b+b^2 = scale^2 = 1, this simplifies some expressions on paper
R = [coeffs(1:2); -coeffs(2) coeffs(1)]/scale; % rotation matrix transforming old->new coordinates (R.' transforms new->old). In particular, R(1,:)*[x;y] is the parametr of the point [x y] when using coeffs_parametric (corresponds to the straight projection (along the axis) onto the parabola)
val1 = coeffs(1)*coeffs(3) + coeffs(2)*coeffs(4);
val2 = coeffs(2)*coeffs(3) - coeffs(1)*coeffs(4);
v = [-coeffs(2) coeffs(1)]/val2;
coeffs_parametric = [v*coeffs(5); (coeffs(1:2)+v*val1)/scale; v*scale2]; % [c0; c1; c2]; follows from writing the paraemtric form in the new coords (y=f(t), x=t) and applying transform new->old
coeffs_new = [scale*coeffs(5) val1 scale^3]/val2; % [c0 c1 c3] coeffs of the parabola y=c2^2+c1*x+c0 in the new coordinates 
end	