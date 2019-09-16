function [f_img m_img state roi] = estimateFM_motion_template(g, b, h, f, m, T, state, varargin)
% Image+mask (=object)-estimation subroutine in the blind esimate-F/estimate-H alternating loop or standalone non-blind image/mask estimation for known PSF. The PSF is a regular motion PSF
%
% g,b - input, background resp., double RGB, same size
% h - psf, same 2-size as g,b, 3rd dimension conatins angles
% f - [] or initial 'f', rgb same 2-size as mask and same channel count as 'g','b'
% m - [] or initial mask (logical)
% One of 'f', 'm', must be spcified to determine size of 'f', 'm'. If more than one is specified, the sizes must match.
% T - template same size as 'f', the penalty is min |f-m*T| (ie teplate is the object at rest, without exaclty specifying the object shape)
% state - strcut containing auxiliary varibles required to resume optimization; simply 'state' as in output

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, @()struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
gamma = IF(isfield(params, 'gamma'), @()params.gamma, 1); % data term weight
alpha = IF(isfield(params, 'alpha'), @()params.alpha, 2^-12); %f tv regularizer weight
alpha_m = IF(isfield(params, 'alpha_m'), @()params.alpha_m, 0); % mask tv regularizer weight
alpha_ml1 = IF(isfield(params, 'alpha_ml1'), @()params.alpha_ml1, 0); % mask l1 regularizer weight
lambda = IF(isfield(params, 'lambda'), @()params.lambda, 0); % template L2 term weight
lambda_m0 = IF(isfield(params, 'lambda_m0'), @()params.lambda_m0, 0); % mask L2 term weight - term lambda*|m-m0|^2 is added to the loss; lambda_m0 can be scalar or 'image' the of 'm' for space-variant weighting |m-m0|^2_w
m0 = IF(isfield(params, 'm0'), @()params.m0, 0); % mask L2 term weight - term lambda*|m-m0|^2 is added to the loss
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 100); % max number of outer iterations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 1e-4); % relative between iterations difference for outer admm loop
cg_maxiter = IF(isfield(params, 'cg_maxiter'), @()params.cg_maxiter, 200); % max number of inner CG iterations ('h' subproblem)
cg_tol = IF(isfield(params, 'cg_tol'), @()params.cg_tol, [1e-4]); % tolerance for relative residual of inner CG iterations ('h' subproblem); can be several values which are used sequentially each time the convergence criterion holds
beta_tv = IF(isfield(params, 'beta_tv'), @()params.beta_tv, 2e1*alpha); % splitting vx/vy=Df due to TV regularizer
beta_tv_m = IF(isfield(params, 'beta_tv_m'), @()params.beta_tv_m, 2e1*alpha_m); % splitting vx_m/vy_m=Dm due to TV regularizer for the mask (equivalent viewpoint - same splitting as vx/y but different alpha/beta param)
beta_ml1 = IF(isfield(params, 'beta_ml1'), @()params.beta_ml1, 0); % splitting m=vml1 due mask L1 regularizer
Lp = IF(isfield(params, 'lp'), @()params.lp, 1); % p-exponent of the "TV" Lp regularizer sum |Df|^p, allowed values are 1/2, 1
Lp_m = IF(isfield(params, 'lp_m'), @()params.lp_m, 1); % p-exponent of the "TV" Lp regularizer sum |Dm|^p for the mask, allowed values are 1/2, 1
beta_fm = IF(isfield(params, 'beta_fm'), @()params.beta_fm, 1e-4); % splitting vf=f adn vm=m due to (f,m) \in C constraint where C is prescribed convex set given by positivity and f-m relation, penalty weight
pyramid_eps = IF(isfield(params, 'pyramid_eps'), @()params.pyramid_eps, 1); % inverse slope of the f<=m/eps constraing for each channel. eps=0 means no constraint (only m in [0,1], f>0), eps=1 means f<=m etc.
fig = IF(isfield(params, 'fig'), @()params.fig, []); % figure handle or id for plotting intermediate results; set to [] to disable plotting
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % 1=display info message in each iteration, 0=don't
do_cost = IF(isfield(params, 'do_cost'), @()params.do_cost, false); % true = calculate cost function value (and related quantities) in each iteration (takes time ~ 1CG iteration); the result it returned in 'cost' struct; false = don't
outer_iter = IF(isfield(params, 'outer_iter'), @()params.outer_iter, 1); % iter in the outer loop, for dbg purposes

if(verbose)
	timer = tic(); % DBG
end
do_cost = do_cost && nargout >= 4; % calculate cost value in each iteration (slower)
rel_tol = rel_tol^2;

% determine sizes of 'f', 'm'
if(~isempty(f))
	fsize = [size(f,1) size(f,2) size(g,3)];
elseif(~isempty(m))
	fsize = [size(m) size(g,3)];
else
	error('One of `f` or `m` must be specified to determine size. Pass trivial value like `zeros` etc.');
end

% crop images to minimal relevant roi (speeds-up the calculation of FTs)
roi = boundingBox(h);
pad = ceil((fsize(1:2)-1)/2); % necessary amount of hmask padding when cropping the input to minimal size required in the FTs
roi = roi + [-pad(1), +pad(1), -pad(2), +pad(2)]; % relevant region in the input image
roi = min(max(roi, [1 1 1 1]), size2(h,[1 1 2 2]));
g = g(roi(1):roi(2), roi(3):roi(4), :);
b = b(roi(1):roi(2), roi(3):roi(4), :);
h = h(roi(1):roi(2), roi(3):roi(4), :);
gsize = size2(g, 1:3);

% unknowns
if(isempty(f))
	f = zeros(prod(fsize(1:2)), fsize(3)); % RGB channels of 'f' a column vectors
else
	f = vec3(f); % RGB channels of 'f' a column vectors
end
if(isempty(m))
	m = ones(prod(fsize(1:2)),1); % intial mask as column vector
else
	m = double(m(:)); % intial mask as column vector
end
if(lambda > 0) T = vec3(T); else T = 0; end % template matching
m0 = double(m0(:)); lambda_m0 = lambda_m0(:); % |m-m0|^2 term and its pixel-wise weighting
F = zeros(gsize); M = zeros(gsize(1:2)); % prealocation for FT calculation (large psfshifted versions of f/m for conv in FT)
[idx_large_f] = psfshift_idx(fsize, gsize); % indices for quick conversion between small 'f' (actual 'f') and 'F' (psf-shifted large version of 'f' for conv in FT)
[idx_large_m] = psfshift_idx(fsize(1:2), gsize(1:2)); % indices for quick conversion between small 'm' (actual 'm') and 'M' (psf-shifted large version of 'm' for conv in FT)

% initial vals
if(nargin >= 8 && ~isempty(state))
	% retrive stored variables
	vx = state.vx; vy = state.vy; ax = state.ax; ay = state.ay; % v=Df splitting due to TV and its assoc. Lagr. mult.
	vx_m = state.vx_m; vy_m = state.vy_m; ax_m = state.ax_m; ay_m = state.ay_m; % v_m=Dm splitting due to TV (m-part) and its assoc. Lagr. mult.
	vf = state.vf; af = state.af; % vf=f splitting due to positivity and f=0 outside mask constraint
	vm = state.vm; am = state.am; % vm=m splitting due to mask between [0,1]
	vml1 = state.vml1; aml1 = state.aml1; % vml1=m splitting due to mask L1 regularizer

	% derivatives
	Dx = state.Dx; Dy = state.Dy; DTD = state.DTD;
else
	vx = zeros(size(f)); vy = zeros(size(f)); ax = zeros(size(f)); ay = zeros(size(f)); % v=Df splitting due to TV and its assoc. Lagr. mult.
	vx_m = zeros(size(m)); vy_m = zeros(size(m)); ax_m = zeros(size(m)); ay_m = zeros(size(m)); % v_m=Dm splitting due to TV (m-part) and its assoc. Lagr. mult.
	vf = zeros(size(f)); af = zeros(size(f)); % vf=f splitting due to positivity and f=0 outside mask constraint
	vm = zeros(size(m)); am = zeros(size(m)); % vm=m splitting due to mask between [0,1]
	vml1 = zeros(size(m)); aml1 = zeros(size(m)); % vml1=m splitting due to mask L1 regularizer
	
	% derivatives
	[Dx Dy] = createDerivatives(fsize); DTD = Dx.'*Dx + Dy.'*Dy;
end

% precompute FT
H = fft2(h); HT = conj(H);

% precompute const RHS for 'f/m' subproblem
rhs_f = real(ifft2(HT.*fft2(g-b))); rhs_f = gamma*reshape(rhs_f(idx_large_f), [], fsize(3));
rhs_m = real(ifft2(HT.*fft2(sum(b.*(g-b),3)))); rhs_m = -gamma*rhs_m(idx_large_m) + lambda_m0.*m0;

% derivatives (FIXME: can be moved to the beginning of admm loop in the relase version and the derivatives calculated before cost function evaluation can be removed)
fdx = Dx*f; fdy = Dy*f;
mdx = Dx*m; mdy = Dy*m;

beta_tv4 = [beta_tv(ones(1,fsize(3))) beta_tv_m];

% DBG - cost function evolution
cost = struct('cost', cell(1,maxiter), 'err', cell(1,maxiter));

% ADMM loop
for iter = 1:maxiter
	f_old = f; m_old = m;
	
	% dual/auxiliary var updates
	% vx/vy minimization (splitting due to TV regularizer)
	if(alpha > 0 && beta_tv > 0)
		val_x = fdx+ax; val_y = fdy+ay;

		[shrink_factor] = lp_proximal_mapping(sqrt(val_x.^2 + val_y.^2), alpha/beta_tv, Lp); % isotropic "TV"
		%[shrink_factor] = lp_proximal_mapping(sqrt(sum(val_x.^2,3) + sum(val_y.^2,3)), alpha/beta_tv, Lp); % isotropic color TV
		% NOTE: anisotropic TV is factor_x = lp_proximal_mapping(abs(val_x), alpha/beta_tv, Lp), the same for _y (ie different shrink factors for val_x and val_y)
				
		vx = val_x.*shrink_factor;
		vy = val_y.*shrink_factor;

		% 'a' step
		ax = ax + fdx - vx;
		ay = ay + fdy - vy;
	end

	% vx_m/vy_m minimization (splitting due to TV regularizer for the mask)
	if(alpha_m > 0 && beta_tv_m > 0)
		val_x = mdx + ax_m; val_y = mdy + ay_m;
		[shrink_factor] = lp_proximal_mapping(sqrt(val_x.^2 + val_y.^2), alpha_m/beta_tv_m, Lp_m); % isotropic TV

		vx_m = val_x.*shrink_factor;
		vy_m = val_y.*shrink_factor;

		% 'a' step
		ax_m = ax_m + mdx - vx_m;
		ay_m = ay_m + mdy - vy_m;
	end
	
	% vf/vm minimization (positivity and and f-m relation so that 'f' is not large for 'm' small)
	if(beta_fm > 0)
		vf = f + af;
		vm = m + am;

		% projection to pyramid-like convex region
		if(pyramid_eps > 0) % (m,f) constrained to convex pyramid-like shape to force f<=const*m where const = 1/eps
			[vm, vf] = project2pyramid(vm, vf, pyramid_eps);
		else % just positivity and m \in [0,1]
			vf(vf < 0) = 0;
			vm(vm < 0) = 0; vm(vm > 1) = 1;
		end

		% lagrange multiplier (dual var) update
		af = af + f - vf;
		am = am + m - vm;
	end

	% vml1 minimization (mask L1 penalty)
	if(beta_ml1 > 0)
		% vml1 - L1 shrinking
		vml1 = m + aml1;
		temp = vml1 < alpha_ml1/beta_ml1;
		vml1(temp) = 0; % also forces positivity
		vml1(~temp) = vml1(~temp) - alpha_ml1/beta_ml1;

		% lagrange multiplier
		aml1 = aml1 + m - vml1;
	end
	
	% f/m step
	rhs1 = rhs_f + beta_tv*(Dx.'*(vx-ax)+Dy.'*(vy-ay)) + beta_fm*(vf-af); % f-part of RHS
	rhs2 = rhs_m + beta_tv_m*(Dx.'*(vx_m-ax_m)+Dy.'*(vy_m-ay_m)) + beta_fm*(vm-am) + beta_ml1*(vml1-aml1); % m-part of RHS
	[fm, cg_iter, cg_residual] = conjgrad(@Ax, [rhs1 rhs2], [f m], cg_tol, cg_maxiter);
	f = fm(:, 1:fsize(3)); m = fm(:, end);

	% NOTE: this is where the convergence tests should be in the release version

	% derivatives (FIXME: move to the beginning of the admm loop in the release version (no cost function eval))
	fdx = Dx*f; fdy = Dy*f;
	mdx = Dx*m; mdy = Dy*m;

	% cost (DBG)
	% FIXME: remove in the release version, calculating data term error costs several uneccessary ffts similar to one CG iteration
	if(do_cost)
		F(idx_large_f) = f; M(idx_large_m) = m;
		cost(iter).err = sum(reshape(real(ifft2(H.*fft2(F)))-b.*real(ifft2(H.*fft2(M)))-(g-b), [], 1).^2);
		cost(iter).cost = gamma/2*cost(iter).err + alpha*sum(sqrt(fdx(:).^2+fdy(:).^2)) + alpha_m*sum(sqrt(mdx(:).^2+mdy(:).^2)) + alpha_ml1*sum(abs(m(:))) + lambda/2*sum(sum((f-m.*T).^2)) + sum(lambda_m0.*(m-m0).^2)/2;
	end

	% FIXME: move to convergence test section
	df = f-f_old; dm = m-m_old;
	rel_diff2_f = (df(:)'*df(:))/(f(:)'*f(:)); % relative l2 norm squared (between-iterations difference)
	rel_diff2_m = (dm(:)'*dm(:))/(m(:)'*m(:)); % relative l2 norm squared (between-iterations difference)

	% DBG
	if(verbose)
		fprintf('FM: iter = %d, cg_iter = %d, cg_res = %.1e, reldiff=(%.1e, %.1e), time=%.2fs\n', iter + outer_iter - 1, cg_iter, cg_residual, sqrt(rel_diff2_f), sqrt(rel_diff2_m), toc(timer));
		timer = tic;
	end
	
	% DBG
	% timer = tic;

	% plot (DBG)
	if(~isempty(fig))
		f_img = ivec3(f, fsize);
		m_img = ivec3(m, fsize);
		figure(fig);
		subplot(121); imshow(f_img); title(sprintf('iter=%d', iter+outer_iter-1)); % title(sprintf('iter=%d, reldiff=%.1e, cgiter=%d, cgresidual=%.1e', iter, sqrt(rel_diff2_f), cg_iter, cg_residual));
		subplot(122); imshow(m_img);% title(sprintf('iter=%d, reldiff=%.1e, cgiter=%d, cgresidual=%.1e', iter, sqrt(rel_diff2_m), cg_iter, cg_residual));
		drawnow;
	end

	% convergence testing
	if(rel_diff2_f < rel_tol && rel_diff2_m < rel_tol)
		break;
	end
end

% reshape images
f_img = ivec3(f, fsize);
m_img = ivec3(m, fsize);

% return state
if(nargout >= 3)
	% return auxiliaries
	state = struct;
	state.vx = vx; state.vy = vy; state.ax = ax; state.ay = ay;
	state.vx_m = vx_m; state.vy_m = vy_m; state.ax_m = ax_m; state.ay_m = ay_m;
	state.vf = vf; state.af = af;
	state.vm = vm; state.am = am;
	state.vml1 = vml1; state.aml1 = aml1;

	% derivatives
	state.Dx = Dx; state.Dy = Dy; state.DTD = DTD;
end

% DBG, cost function value
cost = cost(1:iter);

function y = Ax(x)
xf = x(:, 1:fsize(3)); xm = x(:, end);

F(idx_large_f) = xf; M(idx_large_m) = xm;

% forward conv
HF = H.*fft2(F);
bHM = b.*real(ifft2(H.*fft2(M)));

% transpose conv
yf = real(ifft2(HT.*(HF - fft2(bHM)))); yf = gamma*reshape(yf(idx_large_f),[],fsize(3));
ym = real(ifft2(HT.*fft2(sum(b.*(bHM - real(ifft2(HF))),3)))); ym = gamma*ym(idx_large_m);

% template matching term
if(lambda > 0)
	yf = yf + lambda*(xf-T.*xm);
	ym = ym + lambda*(sum(T.^2,2).*xm - sum(T.*xf,2));
end

% regularizers/identity terms
y = [yf ym + beta_ml1*xm + lambda_m0.*xm] + (DTD*x).*beta_tv4 + x*beta_fm;
end
end

% helper function for vector version of soft thresholding for l1 (lp) minimization
function [shrink_factor] = lp_proximal_mapping(val_norm, amount, p)
shrink_factor = zeros(size(val_norm));
if(p == 1) % soft thresholding
	nz = val_norm > amount;
	shrink_factor(nz) = (val_norm(nz)-amount)./val_norm(nz);
elseif(p == 1/2) % see eg "Computing the proximity operator of the lp norm..., Chen et al, IET Signal processing, 2016"
	nz = val_norm > 3/2*(val_norm).^(2/3);
	shrink_factor(nz) = (2/3*val_norm(nz).*(1+cos(2/3*acos(-3^(3/2)/4*amount*val_norm(nz).^(-3/2)))))./(val_norm(nz));
else
	error('not implemented');
end
end
