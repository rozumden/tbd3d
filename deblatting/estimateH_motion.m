function [H state cost] = estimateH_motion(g, b, f, m, h, hmask, state, varargin)
% Psf-estimation subroutine in the blind esimate-F/estimate-H alternating loop or standalone "non-blind" PSF estimation for known object appearance.
% Simpler version of estimateH_rot, estimates PSF as a classical motion PSF with no rotation
%
% g,b - input, background resp., double RGB, same size (=output size of 'h')
% f - foreground image, same 2-size as 'm'; NOTE: this should actually be f.*m as it appears in the equations, to maintain compatibility with estimateF* and others in the f-h alternating loop (see comment in the beginning of estimateF, eg)
% m - mask (double) of non-transparent pixels in f, same 2-size as 'f'
% h - initial estimate, psf of same 2-size as 'g'. Can be [] as 'none'. Only values in the support of 'hmask' play.
% hmask - logical mask of feasible non-zero points in the translation part of 'h' (same 2-size as 'g' and 'b'). Values outside hmask will be strictly zero - hard constraint. Effectively limits the motion of the object in the image, this restriction does not affect rotation (ie is applied for all angles).
% state - internal; struct containing auxiliary varibles required to resume optimization; simply 'state' as in output; NOTE: if reuse of state is utilized, other factors like 'hmask' or 'roi' (therefore, fsize) or b,g can change, because all auxiliary variables are maintained relative to these variables
%
% cost - struct (1,iter) with cost function calculated in each iteration:
%	.cost - full cost function value excl infty-penalized terms
%	.err - like data term val(^2), but without 1/2 and gamma

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, @()struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
gamma = IF(isfield(params, 'gamma'), @()params.gamma, 1); % data term weight
alpha_h = IF(isfield(params, 'alpha_h'), @()params.alpha_h, 1); % Lp regularizer weight
beta_h = IF(isfield(params, 'beta_h'), @()params.beta_h, 1e3*alpha_h); % splitting for Lp, positivity and sum=1 constraints
% beta_pos = IF(isfield(params, 'beta_pos'), @()params.beta_pos, 0*1e2); % splitting v_pos=h due to positivity constraint; makes sense only for Lp=2, otherwise positivity is enforced by other means without additional splitting
Lp = IF(isfield(params, 'lp'), @()params.lp, 1); % exponent of the Lp regularizer sum |h|^p, allowed values are 0, 1/2, 1, 2
sum1 = IF(isfield(params, 'sum1'), @()params.sum1, true); % force sum(h)=1 constraint (via beta_h). Takes precedence over lp=1 (L1 is ignored, does not make sense)
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 50); % max number of outer iterations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 0); % relative between iterations difference for outer admm loop
cg_maxiter = IF(isfield(params, 'cg_maxiter'), @()params.cg_maxiter, 25); % max number of inner CG iterations ('h' subproblem)
cg_tol = IF(isfield(params, 'cg_tol'), @()params.cg_tol, 1e-5); % tolerance for relative residual of inner CG iterations ('h' subproblem); can be several values which are used sequentially each time the convergence criterion holds
fig = IF(isfield(params, 'fig'), @()params.fig, []); % figure handle or id for plotting intermediate results; set to [] to disable plotting
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % 1=display info message in each iteration, 0=don't
do_cost = IF(isfield(params, 'do_cost'), @()params.do_cost, false); % true = calculate cost function value (and related quantities) in each iteration (takes time ~ 1CG iteration); the result it returned in 'cost' struct; false = don't
outer_iter = IF(isfield(params, 'outer_iter'), @()params.outer_iter, 1); % iter in the outer loop, for dbg purposes

if(verbose)
	timer = tic(); % DBG
end

rel_tol = rel_tol^2;
beta_l2 = 0; % default for lp!=2
if(Lp == 2)
	beta_l2 = alpha_h;
	alpha_h = alpha_h/2; % only appears in cost function calculation
end
do_cost = do_cost && nargout >= 3; % calculate cost value in each iteration (slower)

% crop images to effective roi
fsize = size(f);
hsize = size2(g);

% NOTE: image cropping experimentally removed (investigating influence on result)
pad = ceil((fsize(1:2)-1)/2); % necessary amount of hmask padding when cropping the input to minimal size required in the FTs
if(isempty(hmask)) % proceed as though hmask were all 1s (except for boundary (to prevent wraparound errors))
	hmask = true(hsize);
	hmask_orig = hmask;
	roi = [1 hsize(1) 1 hsize(2)];
else % take hmask as is
	hmask_orig = hmask;
	roi = boundingBox(hmask);
	roi = roi + [-pad(1), +pad(1), -pad(2), +pad(2)];
	roi = min(max(roi, [1 1 1 1]), hsize([1 1 2 2])); % if this line changes 'roi' then wraparound errors may occur (psf is alowed nonzero too close to the boundary)
	% ~~NOTE: ^^ this may later result in error. If the only nz points in hmask are close to the boundary, previous line 'clips' the required neighborhood and subsequent psfshift throws error because fg is then larger then bg (nonsensical) - calling est_H with such input is illegal because relying on wraparound artifacts is unphysical and therefore not correctly supported
	temp = max(fsize(1:2) - [roi(2)-roi(1) roi(4)-roi(3)] - 1, 0); % size deficit
	if(any(temp)) % increase roi size to accomodate 'f' (size can be < due to previous clipping)
		temp2 = min(roi([1 3])-1, temp); % amount of space in top-left direction
		roi([1 3]) = roi([1 3]) - temp2; temp = temp - temp2; % remaining deficit
		roi([2 4]) = roi([2 4]) + temp; % if this causes invalid roi, 'h' is smaller than 'f' and nothing can be done
	end
	g = g(roi(1):roi(2), roi(3):roi(4), :);
	b = b(roi(1):roi(2), roi(3):roi(4), :);
	hmask = hmask(roi(1):roi(2), roi(3):roi(4));
	hsize = size(hmask);
end

% initial values
% note: h contains all unknown (potentially nonzero) pixels in 1 column
if(~isempty(h))
	h = h(roi(1):roi(2), roi(3):roi(4)); h = h(hmask); % recalculate to the roi, if any
else % no initial 'h'
	h = zeros(nnz(hmask),1);
end
if(~isempty(state))
	% retreive stored variables
	v_lp = state.v_lp; % used for lp, positivity and sum=1
	a_lp = state.a_lp;
	% v_pos = state.v_pos;
	% a_pos = state.a_pos;
	% Fgb = state.Fgb; % in the RHS of h-problem, const during subsequent calls and can be stored in the 'state'
	% Fbgb = state.Fbgb;
else
	v_lp = 0; a_lp = 0;
	% v_pos = 0; a_pos = 0;
	% Fgb = fft2(g-b); % in the RHS of h-problem, const during subsequent calls and can be stored in the 'state'
	% Fbgb = fft2(b.*(g-b));
end

% precompute RHS for the 'h' subproblem
F = fft2(psfshift(f, hsize));
M = fft2(psfshift(m, hsize));
Fgb = fft2(g-b);
Fbgb = fft2(b.*(g-b));
rhs_const = sum(real(ifft2(conj(F).*Fgb-conj(M).*Fbgb)),3);
rhs_const = gamma*rhs_const(hmask);

% preallocation for CG iterations
H = zeros(hsize);

% DBG - cost function evolution
cost = struct('cost', cell(1,maxiter), 'err', cell(1,maxiter));

% ADMM loop
for iter = 1:maxiter
	h_old = h;

	% 'v_lp' minimization
	if(beta_h > 0) % either lp, positivity or sum=1 constraint is applied
		v_lp = h + a_lp;

		if(alpha_h > 0 || ~sum1) % lp or positivity
			% Lp
			if(Lp == 1 && ~sum1) % sum=1 constriant implies L1 minimization is senseless
				temp = v_lp < alpha_h/beta_h;
				v_lp(temp) = 0; % also forces positivity
				v_lp(~temp) = v_lp(~temp) - alpha_h/beta_h;
			elseif(Lp == .5)
				% see eg "Computing the proximity operator of the lp norm..., Chen et al, IET Signal processing, 2016"
				temp = v_lp <= 3/2*(alpha_h/beta_h)^(2/3);
				v_lp(temp) = 0; % also forces positivity
				v_lp(~temp) = 2/3*v_lp(~temp).*(1+cos(2/3*acos(-3^(3/2)/4*alpha_h/beta_h*v_lp(~temp).^(-3/2))));
			elseif(Lp == 0)
				v_lp(v_lp <= sqrt(2*alpha_h/beta_h)) = 0; % also forces positivity
			end
		end

		if(sum1) % sum=1 constraint
			v_lp = proj2simplex(v_lp);
		end

		% 'a' step
		a_lp = a_lp + h - v_lp;
	end
	
	% % 'v_pos' minimization
	% if(beta_pos > 0)
	% 	v_pos = h + a_pos;
	% 	v_pos(v_pos < 0) = 0; % infinit penalty for negative 'h' values
	% 	v_pos(v_pos > 1) = 1; % infinit penalty for 'h' values over 1 
		
	% 	% 'a' step
	% 	a_pos = a_pos + h - v_pos;
	% end

	% 'h' minimization
	% rhs = rhs_const + beta_pos*(v_pos-a_pos) + beta_h*(v_lp-a_lp); % last term is zero for lp=2
	rhs = rhs_const + beta_h*(v_lp-a_lp);

	% CG solution
	[h, cg_iter, cg_residual] = conjgrad(@Ax, rhs, h, cg_tol, cg_maxiter);

	% NOTE: this is where the convergence tests should be in the release version
	
	% cost (DBG)
	if(do_cost)
		% note: calculating data term error costs several uneccessary ffts similar to one CG iteration
		H(hmask) = h; % unpack 'h' to the 3d psf
		FH = fft2(H);
		% apply forward conv (->RGB image, summation over angles)
		Fh = F.*FH;
		BMh = b.*real(ifft2(M.*FH));
		cost(iter).err = sum(reshape(real(ifft2(Fh))-BMh-(g-b), [], 1).^2);
		cost(iter).cost = gamma/2*cost(iter).err + alpha_h*sum(reshape(abs(h).^Lp, [], 1));
	end
	
	% convergence measure
	d = h-h_old;
	rel_diff2 = (d(:)'*d(:))/(h(:)'*h(:)); % relative l2 norm squared (between-iterations difference)
	
	% info message
	if(verbose)
		fprintf('H: iter=%d, cost=(%.1e/%.1e), cg_iter=%d, cg_res=%.1e, reldiff=%.1e, time=%.2fs\n', iter, cost(iter).err, cost(iter).cost, cg_iter, cg_residual, sqrt(rel_diff2), toc(timer));
		timer = tic();
	end

	% plot (DBG)
	if(~isempty(fig))
		h_temp = zeros(size(hmask_orig)); h_temp(hmask_orig) = h;
		figure(fig);
		imshow(h_temp/max(h)); title(sprintf('iter=%d, sum=%.2f, sum(abs)=%.2f', iter+outer_iter-1, sum(h), sum(abs(h))));
		drawnow;
	end

	% convergence testing
	if(rel_diff2 < rel_tol)
		break;
	end
end

% return as full-sized psf, pad back what was cropped by roi
H = zeros(size(hmask_orig));
H(hmask_orig) = h;

% return state
if(nargout >= 2)
	state = struct;
	state.v_lp = v_lp;
	state.a_lp = a_lp;
	% state.v_pos = v_pos;
	% state.a_pos = a_pos;
	% state.Fgb = Fgb;
	% state.Fbgb = Fbgb;
end

% DBG, cost function value
cost = cost(1:iter);

function res = Ax(h) % inner function for the CG iterations
H(hmask) = h; % unpack 'h' to the 3d psf
FH = fft2(H);

% apply forward conv (->RGB image, summation over angles)
Fh = F.*FH;
BMh = b.*real(ifft2(M.*FH));
Fh_BMh = Fh - fft2(BMh); % FT of Fh - BMh
% transposed conv
res = sum(real(ifft2(conj(F).*Fh_BMh - conj(M).*fft2(b.*real(ifft2(Fh_BMh))))),3);

% identity terms (regularization)
% res = gamma*res(hmask) + (beta_h+beta_pos)*h; % identity terms, with v_lp splitting (naturally includes lp=2 with beta_h:=alpha_h)
res = gamma*res(hmask) + (beta_h+beta_l2)*h; % identity terms
end
end

function x = proj2simplex(y)
% euclidean projection of y (arbitrarily shaped but treated as a single vector) to a simplex defined as x>=0 and sum(x(:)) = 1
% based on "Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application"; Weiran Wang et al; 2013 (arXiv:1309.1541)
x = sort(y(:),'descend');
temp = (cumsum(x)-1)./(1:numel(x)).';
x = reshape(max(y(:) - temp(find(x > temp, 1, 'last')), 0), size(y));
end
