function [coeffs control_pts dist] = fitPwPolynomial(pts, control_pts, num_curves, w, maxiter, tol, sampling_rate, epsilon)
% Fits a continuous parametric polynomial curve of specified order to the set of points pts. The curve is broken into num_curves segments, each being a polynomial of the given order, but togenther conencted only C^0.
% Different from fitPoly in 2 ways: 1/ allows breaking the polynomial into multiple segments, 2/ performs minimization wrt control points and not coefficients (probably doesn't matter). This is more mature implementation than fitPoly.
% The polynomial order is given by the initialization, ie number of control points.
% The composite curve is given as several (num_curves) mappings t:[0,1]->R^d of the form point(t)=c(0)+c(1)*t+c(2)*t^2+c(3)*t^3+... where c(i) = coeffs(i,:). (ie, no 1/k! factors are included). It is guaranteed that curve_k(1)=curve_{k+1}(0), ie curves are connected.
%
% pts - input points to be fitted, Nxd set of N points with 'd' coordinates
% control_pts - initialization; Mxd array of polynomial control points, M is the polynomial order. Control points are attained values at equidistant (in terms of parametrization) times for each curve, last control points of one curve and first of the next curve are shared.
%  - the number of control points must be N*num_curves+1 where N is the polynomial order.
% num_curves - numebr of curve segments (polynomial curves), default=1
% w - optional weights of points vector, same length as number of input pts
% maxiter - max number of iterations
% tol - stopping tolerance; between iteration change - mean (or max, see the code) displacement of any control point between iterations
% sampling_rate - number of discretized steps in the fitted curve (each segment) segment sampling
% epsilon - irls epsilon (roughly, loss is qudratic for |val| < eps, linear otherwise (or similar, see code)); set to inf for fully quadratic loss (=disable irls)
%
%
% coeffs - cell array (1,num_curves), each containing polynomial coefficients of the particular curve. Indexed as (power increasing top->bottom,coordinate) = similar format as control_pts and same size for each particular curve
% control_pts - all control pts of the composite curve together
% dist - distances of input pts to the curve

% default args
if(nargin < 3) num_curves = 1; end
if(nargin < 4) w = ones(size(pts,1),1); end
if(nargin < 5) maxiter = 20*num_curves; end
if(nargin < 6) tol = 1e-3; end
if(nargin < 7) sampling_rate = 100; end
if(nargin < 8) epsilon = .5; end

% polynomial order N: N+1 points are for each curve, so N*num_curves+1 is the total number of points, since first/last of each curve are shared
N = (size(control_pts,1)-1)/num_curves;

% hard-coded params
% approx_len = sum(sqrt(sum((control_pts(2:end,:)-control_pts(1:end-1,:)).^2,2)))/num_curves; % mean length of each curve (lower bound)
% sampling_rate = 2*ceil(approx_len); % number of discretized steps in the fitted curve (each segment) segment sampling
tie_ends = true; % restrict curve ends to particular input points (restricts potentialy infinite growth of the curve). 0=no restriction, 1= means that first and last point (as projected to the curve) will be measured against curve beginning/endresp.. Note that this does not mean that the curve endpoints are fixed.
control_t = linspace(0,1,N+1); % parametr values (times) corresponding to given control pts. Must not contain (close to) duplicate values, otherwise the corresponding Vandermonde matrix is (close to) singular. Common for all segments and each segment is parametrized [0,1]->, so spcify these accordingly.

% precompute sampling matrix
V = (control_t(:).^(0:N)); % Vandermonde matrix relating polynomial coefficients and control point in the given-degree polynomial: control_pts = V*polynomial_coefficients. Common for all curves as long as control_t and of course orders are same.
tv = ((linspace(0,1,sampling_rate+1).').^(0:N))/V; % precomputed coeffs that appear in the gradient of the approximated problem and in curve sampling, indexed as (time_idx, power)

% prealocations and init
num_pts = size(tv,1); % number of sampled points on each curve segment
curve_pts = zeros(num_curves*num_pts,size(pts,2)); % curve pts of all segments; 1:num_pts is first curve, num_pts+1:2*num_pts is second curve etc
lhs = zeros(size(control_pts,1),size(control_pts,1)); % LHS matrix of the linear problem in each iter; contains subproblems for all curve segments with the eqs for shared points summed
rhs = zeros(size(control_pts)); % RHS for the main linear system. Again, eqs for all control points incl. sharing
w_irls = ones(size(pts,1),1);

% main loop
for iter=1:maxiter
	old_pts = control_pts;
	
	% trace curve
	for i=1:num_curves
		curve_pts((i-1)*num_pts+1:i*num_pts,:) = tv*control_pts((i-1)*N+1:i*N+1,:); % indexed as (point_idx,coordinate) where point idx is actuall curve_idx*point_idx; N is the offset in control points between consecutive curves (same as polynomial order, number of control pts in a curve-1)
	end

	% find points matching (=for each input point, find matching curve point)
	dist2 = sum((reshape(pts,size(pts,1),1,size(pts,2))-reshape(curve_pts,1,size(curve_pts,1), size(curve_pts,2))).^2,3); % (squared) distance "input_point-curve_point", all mutual combinations; indexed as (input_pts,curve_pts)
	[dist2,idx] = min(dist2,[],2); % closest curve point to each input point; idx needs to be recalculated to curve_idx, time_idx
	
	% modification of loss to measure particular points to curve ends (restricts curve growth to infty)
	if(tie_ends)
		[~,i1] = min(idx); [~,i2] = max(idx);
		idx([i1 i2]) = [1 num_curves*num_pts]; % tie first point to beginning of first curve, last point to end of last curve
	end

	% DBG
% 	figure(1);
% 	plot(pts(:,2), pts(:,1), 'k+'); title(sprintf('iter=%d, loss=%.2e', iter, mean(sqrt(dist2(:)))));
% 	hold on;
% 	plot([pts(:,2).'; curve_pts(idx,2).'], [pts(:,1).'; curve_pts(idx,1).'], '-ro'); % matching pts
% 	plot(curve_pts(:,2), curve_pts(:,1), '-b');
% 	hold off;
% 	axis equal;
% 	drawnow;
% 	pause;

	% recalculate idx to curve/point
	pt_idx = mod(idx-1,num_pts)+1; % index into particular point in the curve
	curve_idx = (idx-pt_idx)/num_pts+1; % 1-based curve index, where the given input pt is projected

	% l1 reweighting
	if(isfinite(epsilon))
		% NOTE: dont' know which is better or 'more correct'
		%w_irls = max(sqrt(dist2), epsilon);
		w_irls = sqrt(dist2+epsilon^2);
	end

	% construct gradient (lhs and rhs of the corresponding linear system "grad wrt coeffs = 0")
	lhs(:) = 0; rhs(:) = 0; % reset system, build equations for each curve (eqs for connecting points are shared(summed))
	for i=1:num_curves
		mask = curve_idx == i;
		span = (i-1)*N+1:i*N+1; % which indices in the linear system are affected by the processed curve
		v = tv(pt_idx(mask),:); vt = v.'; % temp, selected (active) samples of the sampling matrix tv
		lhs(span, span) = lhs(span, span) + vt*(v.*w(mask)./w_irls(mask)); % add equations for particular curve
		rhs(span,:) = rhs(span,:) + vt*(pts(mask,:).*w(mask)./w_irls(mask));
	end
	control_pts = lhs\rhs; % solve corresponding system

	% stopping criterion - max or average displacement of any control pt
	diff = max(sum((control_pts-old_pts).^2,2)); % maximum between iteration displacement of any control point (squared)
	if(diff < tol^2)
		break;
	end
end

% return polynomial coeffs of each segment separately
coeffs = cell(1,num_curves);
for i=1:num_curves
	coeffs{i} = V\control_pts((i-1)*N+1:i*N+1, :); % V^-1 transforms control pts to coefficients
end

% return distances if required
if(nargout >= 3)
	% trace curve
	for i=1:num_curves
		curve_pts((i-1)*num_pts+1:i*num_pts,:) = tv*control_pts((i-1)*N+1:i*N+1,:); % indexed as (point_idx,coordinate) where point idx is actuall curve_idx*point_idx; N is the offset in control points between consecutive curves (same as polynomial order, number of control pts in a curve-1)
	end

	% find points matching (=for each input point, find matching curve point)
	dist = sum((reshape(pts,size(pts,1),1,size(pts,2))-reshape(curve_pts,1,size(curve_pts,1), size(curve_pts,2))).^2,3); % (squared) distance "input_point-curve_point", all mutual combinations; indexed as (input_pts,curve_pts)
	[dist,idx] = min(dist,[],2); % closest curve point to each input point; idx needs to be recalculated to curve_idx, time_idx
	dist = sqrt(dist);
end
end