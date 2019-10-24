function [f_img m_img state cost] = estimateFM_motion_pw(g, b, h, f, m, T, state, varargin)
% Image+mask (=object)-estimation subroutine in the blind esimate-F/estimate-H alternating loop or standalone non-blind image/mask estimation for known PSF.
% The model is essentialy "multiple objects+multiple PSFs", though mainly intended as piecewise-const apearance of one object and one motion. The psf is stacked PSFs in 3-dim (ie sum(h,3) is the whole PSF) and F and M are corresponding image/mask for individual segments (optionally common mask for all objects)
%
% g,b - input, background resp., double RGB, same size
% h - psf, same 2-size as g,b, 3rd dimension is different parts of the piecewise const-appearance motion or different psfs of correpsonding objects
% f - [] or initial 'f', rgb same 2-size as mask and same channel count as 'g','b'; different appearances stacked in 4th dim as HWCN; must be spcified to determine size of 'f', 'm'
% m - [] or initial mask; ; different appearances stacked in 3rd or 4th dim; NOTE: when initial m of size HW11 is passed, then different implementation is used - single common 'm' for all objects
% One of 'f', 'm', must be spcified to determine size of 'f', 'm'. If more than one is specified, the sizes must match.
% T - template for 'F', constriant is lambda*|F_i - T_i|^2 (ie - N templates can be specified for each F_i); if lambda != 0 and and F!= [] and T=[] then T is initializaed as T_i:=F_i; the code contains commented-out version with the original formulation lambda*|F_i-M_i*T_i| - requires modifying Ax and rhs_f
% state - strcut containing auxiliary varibles required to resume optimization; simply 'state' as in output

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, @()struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
gamma = IF(isfield(params, 'gamma'), @()params.gamma, 1); % data term weight
alpha_f = IF(isfield(params, 'alpha_f'), @()params.alpha_f, 2^-12); %f tv regularizer weight
alpha_m = IF(isfield(params, 'alpha_m'), @()params.alpha_m, 2^-12); % mask tv regularizer weight
alpha_cross_f = IF(isfield(params, 'alpha_cross_f'), @()params.alpha_cross_f, 0); % cross-image (in 3-dim) image tv regularizer weight
alpha_cross_m = IF(isfield(params, 'alpha_cross_m'), @()params.alpha_cross_m, 0); % cross-image (in 3-dim) mask tv regularizer weight (ignored if single_m)
lambda = IF(isfield(params, 'lambda'), @()params.lambda, 0); % template L2 term weight (scalar or matrix same size2 as template for spacially-variant weighting)
lambda_R = IF(isfield(params, 'lambda_R'), @()params.lambda_R, 0); % mask rotation symmetry weight term, lambda_R*|R*m-m|^2 where R is approx rotational averaging (applied to mask only)
lambda_m0 = IF(isfield(params, 'lambda_m0'), @()params.lambda_m0, 0); % mask L2 term weight - term lambda*|m-m0|^2 is added to the loss; lambda_m0 can be scalar or 'image' the of 'm' for space-variant weighting |m-m0|^2_w
m0 = IF(isfield(params, 'm0'), @()params.m0, 0); % mask L2 term weight - term lambda*|m-m0|^2 is added to the loss
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 50); % max number of outer iterations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 0); % relative between iterations difference for outer admm loop
cg_maxiter = IF(isfield(params, 'cg_maxiter'), @()params.cg_maxiter, 25); % max number of inner CG iterations ('fm' subproblem)
cg_tol = IF(isfield(params, 'cg_tol'), @()params.cg_tol, 1e-5); % tolerance for relative residual of inner CG iterations ('fm' subproblem); can be several values which are used sequentially each time the convergence criterion holds
beta_f = IF(isfield(params, 'beta_f'), @()params.beta_f, 10*alpha_f); % splitting vx/vy=Df due to TV regularizer
beta_m = IF(isfield(params, 'beta_m'), @()params.beta_m, 10*alpha_m); % splitting vx_m/vy_m=Dm due to TV regularizer for the mask (equivalent viewpoint - same splitting as vx/y but different alpha_f/beta param)
beta_cross_f = IF(isfield(params, 'beta_cross_f'), @()params.beta_cross_f, 10*alpha_cross_f); % splitting vc=D_cross*f due to cross-image TV regularizer
beta_cross_m = IF(isfield(params, 'beta_cross_m'), @()params.beta_cross_m, 10*alpha_cross_m); % splitting vc_m=D_cross*m due to cross-mask TV regularizer (not used is single_m)
Lp = IF(isfield(params, 'lp'), @()params.lp, 1); % p-exponent of the "TV" Lp regularizer sum |Df|^p, allowed values are 1/2, 1
Lp_m = IF(isfield(params, 'lp_m'), @()params.lp_m, 1); % p-exponent of the "TV" Lp regularizer sum |Dm|^p for the mask, allowed values are 1/2, 1
beta_fm = IF(isfield(params, 'beta_fm'), @()params.beta_fm, 1e-3); % splitting vf=f adn vm=m due to (f,m) \in C constraint where C is prescribed convex set given by positivity and f-m relation, penalty weight
pyramid_eps = IF(isfield(params, 'pyramid_eps'), @()params.pyramid_eps, 1); % inverse slope of the f<=m/eps constraing for each channel. eps=0 means no constraint (only m in [0,1], f>0), eps=1 means f<=m etc.
fig = IF(isfield(params, 'fig'), @()params.fig, []); % figure handle or id for plotting intermediate results; set to [] to disable plotting
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % 1=display info message in each iteration, 0=don't
do_cost = IF(isfield(params, 'do_cost'), @()params.do_cost, false); % FIXME: not implemented; true = calculate cost function value (and related quantities) in each iteration (takes time ~ 1CG iteration); the result it returned in 'cost' struct; false = don't
outer_iter = IF(isfield(params, 'outer_iter'), @()params.outer_iter, 1); % iter in the outer loop, for dbg purposes

if(verbose)
	timer = tic(); % DBG
end
rel_tol = rel_tol^2;

% determine sizes of 'f', 'm'
if(~isempty(f))
	fsize = [size(f,1) size(f,2) size(g,3) size(h,3)];
elseif(~isempty(m))
	fsize = [size2(m) size(g,3) size(h,3)];
else
	error('One of `f` or `m` must be specified to determine size. Pass trivial value like `zeros` etc.');
end

% crop images to minimal relevant roi (speeds-up the calculation of FTs)
roi = boundingBox(any(h,3));
pad = ceil((fsize(1:2)-1)/2); % necessary amount of hmask padding when cropping the input to minimal size required in the FTs
roi = roi + [-pad(1), +pad(1), -pad(2), +pad(2)]; % relevant region in the input image
roi = min(max(roi, [1 1 1 1]), size2(h,[1 1 2 2]));
temp = max(fsize(1:2) - [roi(2)-roi(1) roi(4)-roi(3)] - 1, 0); % size deficit
if(any(temp)) % increase roi size to accomodate 'f' (size can be < due to previous clipping)
	temp2 = min(roi([1 3])-1, temp); % amount of space in top-left direction
	roi([1 3]) = roi([1 3]) - temp2; temp = temp - temp2; % remaining deficit
	roi([2 4]) = roi([2 4]) + temp; % if this causes invalid roi, 'h' is smaller than 'f' and nothing can be done
end
g = g(roi(1):roi(2), roi(3):roi(4), :);
b = b(roi(1):roi(2), roi(3):roi(4), :);
h = h(roi(1):roi(2), roi(3):roi(4), :);

% unknowns
if(isempty(f))
	f = zeros(prod(fsize(1:2)), fsize(3)*fsize(4)); % RGB channels of 'f' a column vectors, different objects(apperences) stacked in 2-dim as [RGB1 RGB2 ...]
else
	f = vec3(f); % RGB channels of 'f' a column vectors, different objects(apperences) stacked in 2-dim as [RGB1 RGB2 ...]
end
if(isempty(m))
	m = ones(prod(fsize(1:2)), fsize(4)); % intial mask as column vector, different objects(apperences) stacked in 2-dim
else
	m = vec3(m); % intial mask as column vector
end
single_m = size(m,2) == 1; % single mask is const for all object appearances
lambda = lambda(:);
if(any(lambda > 0))
	if(~isempty(T)) T = vec3(T); else T=f; end % default
else
	lambda = 0; T = 0;
end
m0 = double(m0(:)); lambda_m0 = lambda_m0(:); % |m-m0|^2 term and its pixel-wise weighting
gsize = size2(g, 1:3);
F = zeros([gsize fsize(4)]); M = zeros([gsize(1:2) 1 size(m,2)]); % prealocation for FT calculation (large psfshifted versions of f/m for conv in FT), incl number of objects in 4th dim
[idx_large_f] = psfshift_idx(fsize, [gsize fsize(4)]); % indices for quick conversion between small 'f' (actual 'f') and 'F' (psf-shifted large version of 'f' for conv in FT)
[idx_large_m] = psfshift_idx([fsize(1:2) 1 size(m,2)], [gsize(1:2) 1 size(m,2)]); % indices for quick conversion between small 'm' (actual 'm') and 'M' (psf-shifted large version of 'm' for conv in FT)
h = reshape(h, [gsize(1:2) 1 fsize(4)]); % same format as F,M in fft2 - different PSFs in 4th dim

% rescale tv regularizations
alpha_f = alpha_f/fsize(4); beta_f = beta_f/fsize(4); lambda = lambda/fsize(4);
if(fsize(4) > 1) alpha_cross_f = alpha_cross_f/(fsize(4)-1); beta_cross_f = beta_cross_f/(fsize(4)-1); end
if(~single_m)
	alpha_m = alpha_m/fsize(4); beta_m = beta_m/fsize(4);
	alpha_cross_m = alpha_cross_m/(fsize(4)-1); beta_cross_m = beta_cross_m/(fsize(4)-1);
	lambda_R = lambda_R/fsize(4);
	lambda_m0 = lambda_m0/fsize(4);
else
	% disable cross-m derivative
	beta_cross_m = 0;
end

% initial vals
if(~isempty(state))
	% retrive stored variables
	vx = state.vx; vy = state.vy; ax = state.ax; ay = state.ay; % v=Df splitting due to TV and its assoc. Lagr. mult.
	vx_m = state.vx_m; vy_m = state.vy_m; ax_m = state.ax_m; ay_m = state.ay_m; % v_m=Dm splitting due to TV (m-part) and its assoc. Lagr. mult.
	vc = state.vc; ac = state.ac; vc_m = state.vc_m; ac_m = state.ac_m; % vc=D_cross*f splitting due to cross-image/mask TV and its assoc. Lagr. mult. (TV works along the 3rd dim of 'f' or 'm', across the different apeparances/masks of the object)
	vf = state.vf; af = state.af; % vf=f splitting due to positivity and f=0 outside mask constraint
	vm = state.vm; am = state.am; % vm=m splitting due to mask between [0,1]

	% derivatives
	Dx = state.Dx; Dy = state.Dy; DTD = state.DTD;

	% mask rotation symmetry
	Rn = state.Rn;
else
	% derivatives (cross-img derivates are calculated manually)
	[Dx Dy] = createDerivatives0(fsize); DTD = Dx.'*Dx + Dy.'*Dy;

	vx = zeros(size(Dx,1),size(f,2)); vy = zeros(size(Dy,1),size(f,2)); ax = 0; ay = 0; % v=Df splitting due to TV and its assoc. Lagr. mult.
	vx_m = zeros(size(Dx,1),size(m,2)); vy_m = zeros(size(Dy,1),size(m,2)); ax_m = 0; ay_m = 0; % v_m=Dm splitting due to TV (m-part) and its assoc. Lagr. mult.
	% vc = zeros([fsize(1)*fsize(2) fsize(3)*(fsize(4)-1)]); ac = zeros(size(vc)); % vc=D_cross*f splitting due to cross-image TV and its assoc. Lagr. mult. (TV works along the 3rd dim of 'f', across the different apeparances of the object)
	% vc_m = zeros([fsize(1)*fsize(2) fsize(4)-1]); ac_m = zeros(size(vc_m)); % vc_m=D_cross*m splitting due to cross-mask TV and its assoc. Lagr. mult. (TV works along the 3rd dim of 'm', across the different masks of the object)
	vc = 0; ac = 0; % vc=D_cross*f splitting due to cross-image TV and its assoc. Lagr. mult. (TV works along the 3rd dim of 'f', across the different apeparances of the object)
	vc_m = 0; ac_m = 0; % vc_m=D_cross*m splitting due to cross-mask TV and its assoc. Lagr. mult. (TV works along the 3rd dim of 'm', across the different masks of the object)
	vf = 0; af = 0; % vf=f splitting due to positivity and f=0 outside mask constraint
	vm = 0; am = 0; % vm=m splitting due to mask between [0,1]

	% mask rotation symmetry
	if(lambda_R)
		Rn = createRnMatrix(fsize(1:2));
		Rn = Rn.'*Rn-Rn.'-Rn+speye(size(Rn)); % matrix that actually apears in the quad term
	else
		Rn = 0;
	end
end

% cross-img derivatives and single image (mask) cases
if(fsize(4) == 1 || beta_cross_f == 0) % single image - disable cross-img derivatives (would cause errors)
	crossDf = @(x)0; crossDf_T = @(x)0; crossDf_DTD = @(x)0;
else
	crossDf = @(x)crossD(x,fsize(3)); crossDf_T = @(x)crossD_T(x,fsize(3)); crossDf_DTD = @(x)crossD_DTD(x,fsize(3));
end
if(single_m || beta_cross_m == 0)
	crossDm = @(x)0; crossDm_T = @(x)0; crossDm_DTD = @(x)0;
else
	crossDm = @(x)crossD(x,1); crossDm_T = @(x)crossD_T(x,1); crossDm_DTD = @(x)crossD_DTD(x,1);
end

% precompute FT
H = fft2(h); HT = conj(H);

% precompute const RHS for 'f/m' subproblem
rhs_f = real(ifft2(HT.*fft2(g-b))); rhs_f = gamma*reshape(rhs_f(idx_large_f), fsize(1)*fsize(2),fsize(3)*fsize(4));
rhs_f = rhs_f + lambda.*T; % template matching term with lambda*|F_i-T_i| formulation (remove for masked formulation)
if(single_m)
	rhs_m = real(ifft2(sum(HT,4).*fft2(sum(b.*(g-b),3))));
else
	rhs_m = real(ifft2(HT.*fft2(sum(b.*(g-b),3))));
end
rhs_m = -gamma*reshape(rhs_m(idx_large_m), fsize(1)*fsize(2),[]) + lambda_m0.*m0;

% derivatives (NOTE: can be moved to the beginning of admm loop in the relase version and the derivatives calculated before cost function evaluation can be removed)
fdx = Dx*f; fdy = Dy*f;
mdx = Dx*m; mdy = Dy*m;
fdc = crossDf(f); mdc = crossDm(m);
beta_tv4 = [beta_f(ones(1,fsize(3)*fsize(4))) beta_m(ones(1,size(m,2)))];

% DBG - cost function evolution
if(do_cost)
	cost = struct('cost', cell(1,maxiter), 'err', cell(1,maxiter));
end

% ADMM loop
for iter = 1:maxiter
	f_old = f; m_old = m;
	
	% dual/auxiliary var updates
	% vx/vy minimization (splitting due to TV regularizer)
	if(beta_f > 0)
		val_x = fdx+ax; val_y = fdy+ay;

		[shrink_factor] = lp_proximal_mapping(sqrt(val_x.^2 + val_y.^2), alpha_f/beta_f, Lp); % isotropic "TV"
		%[shrink_factor] = lp_proximal_mapping(sqrt(sum(val_x.^2,3) + sum(val_y.^2,3)), alpha_f/beta_f, Lp); % isotropic color TV
		% NOTE: anisotropic TV is factor_x = lp_proximal_mapping(abs(val_x), alpha_f/beta_f, Lp), the same for _y (ie different shrink factors for val_x and val_y)
				
		vx = val_x.*shrink_factor;
		vy = val_y.*shrink_factor;

		% 'a' step
		ax = ax + fdx - vx;
		ay = ay + fdy - vy;
	end

	% vx_m/vy_m minimization (splitting due to TV regularizer for the mask)
	if(beta_m > 0)
		val_x = mdx + ax_m; val_y = mdy + ay_m;
		[shrink_factor] = lp_proximal_mapping(sqrt(val_x.^2 + val_y.^2), alpha_m/beta_m, Lp_m); % isotropic TV

		vx_m = val_x.*shrink_factor;
		vy_m = val_y.*shrink_factor;

		% vx_m = sign(vx_m).*max(abs(val_x)-alpha_m/beta_m, 0);
		% vy_m = sign(vy_m).*max(abs(val_y)-alpha_m/beta_m, 0);

		% 'a' step
		ax_m = ax_m + mdx - vx_m;
		ay_m = ay_m + mdy - vy_m;
	end

	% cross-image derivative
	if(beta_cross_f > 0)
		val = fdc + ac;

		[shrink_factor] = lp_proximal_mapping(val, alpha_cross_f/beta_cross_f, 1);
		vc = val.*shrink_factor;

		% 'a' step
		ac = ac + fdc - vc;
	end

	% cross-mask derivative
	if(beta_cross_m > 0)
		val = mdc + ac_m;

		[shrink_factor] = lp_proximal_mapping(val, alpha_cross_m/beta_cross_m, 1);
		vc_m = val.*shrink_factor;

		% 'a' step
		ac_m = ac_m + mdc - vc_m;
	end
	
	% vf/vm minimization (positivity and and f-m relation so that 'f' is not large for 'm' small)
	if(beta_fm > 0)
		vf = f + af;
		vm = m + am;

		% projection to pyramid-like convex region
		if(pyramid_eps > 0) % (m,f) constrained to convex pyramid-like shape to force f<=const*m where const = 1/eps
			if(single_m)
				[vm, vf] = project2pyramid(vm, vf, pyramid_eps);
			else
				% note: simple looped version
				for i=1:fsize(4)
					[vm(:,i), vf(:,(i-1)*fsize(3)+1:(i-1)*fsize(3)+3)] = project2pyramid(vm(:,i), vf(:,(i-1)*fsize(3)+1:(i-1)*fsize(3)+3), pyramid_eps);
				end
			end
		else % just positivity and m in [0,1]
			vf(vf < 0) = 0;
			vm(vm < 0) = 0; vm(vm > 1) = 1;
		end

		% lagrange multiplier (dual var) update
		af = af + f - vf;
		am = am + m - vm;
	end
	
	% f/m step
	rhs1 = rhs_f + beta_f*(Dx.'*(vx-ax)+Dy.'*(vy-ay)) + beta_cross_f*crossDf_T(vc-ac) + beta_fm*(vf-af); % f-part of RHS
	rhs2 = rhs_m + beta_m*(Dx.'*(vx_m-ax_m)+Dy.'*(vy_m-ay_m)) + beta_cross_m*crossDm_T(vc_m-ac_m) + beta_fm*(vm-am); % m-part of RHS
	[fm, cg_iter, cg_residual] = conjgrad(@Ax, [rhs1 rhs2], [f m], cg_tol, cg_maxiter);
	f = fm(:, 1:fsize(3)*fsize(4)); m = fm(:, fsize(3)*fsize(4)+1:end);

	% NOTE: this is where the convergence tests should be in the release version

	% derivatives (note: move to the beginning of the admm loop in the release version (no cost function eval))
	fdx = Dx*f; fdy = Dy*f;
	mdx = Dx*m; mdy = Dy*m;
	fdc = crossDf(f); mdc = crossDm(m);


	% cost (DBG)
	% note: remove in the release version, calculating data term error costs several uneccessary ffts similar to one CG iteration
	if(do_cost)
		% FIXME: not updated to _pw version incl cross-img derivatives and lambda terms
		% F(idx_large_f) = f; M(idx_large_m) = m;
		% cost(iter).err = sum(reshape(real(ifft2(H.*fft2(F)))-b.*real(ifft2(H.*fft2(M)))-(g-b), [], 1).^2);
		% cost(iter).cost = gamma/2*cost(iter).err + alpha_f*sum(sqrt(fdx(:).^2+fdy(:).^2)) + alpha_m*sum(sqrt(mdx(:).^2+mdy(:).^2)) + alpha_ml1*sum(abs(m(:)));
	end

	% note: move to convergence test section
	df = f-f_old; dm = m-m_old;
	rel_diff2_f = (df(:)'*df(:))/(f(:)'*f(:)); % relative l2 norm squared (between-iterations difference)
	rel_diff2_m = (dm(:)'*dm(:))/(m(:)'*m(:)); % relative l2 norm squared (between-iterations difference)

	% DBG
	if(verbose)
		fprintf('FM: iter = %d, cg_iter = %d, cg_res = %.1e, reldiff=(%.1e, %.1e), time=%.2fs\n', iter + outer_iter - 1, cg_iter, cg_residual, sqrt(rel_diff2_f), sqrt(rel_diff2_m), toc(timer));
		timer = tic;
	end
	
	% plot (DBG)
	if(~isempty(fig))
		f_img = ivec3(f, fsize);
		m_img = ivec3(m, fsize(1:2));
		figure(fig);
		for i=1:fsize(4)
			subplot(2,fsize(4),i); imshow(f_img(:,:,:,i)); % title(sprintf('iter=%d, reldiff=%.1e, cgiter=%d, cgresidual=%.1e', iter, sqrt(rel_diff2_f), cg_iter, cg_residual));
			if(i == 1) title(sprintf('iter=%d', iter)); end
			subplot(2,fsize(4),fsize(4)+i); imshow(m_img(:,:,min(i, size(m_img,3)))); % title(sprintf('iter=%d, reldiff=%.1e, cgiter=%d, cgresidual=%.1e', iter, sqrt(rel_diff2_m), cg_iter, cg_residual));
		end
		drawnow;
	end

	% convergence testing
	if(rel_diff2_f < rel_tol && rel_diff2_m < rel_tol)
		break;
	end
end

% reshape images
f_img = ivec3(f, fsize);
m_img = ivec3(m, fsize(1:2)); % stacked as 3d, unlike f_img

% return state
if(nargout >= 3)
	% return auxiliaries
	state = struct;
	state.vx = vx; state.vy = vy; state.ax = ax; state.ay = ay;
	state.vx_m = vx_m; state.vy_m = vy_m; state.ax_m = ax_m; state.ay_m = ay_m;
	state.vc = vc; state.vc_m = vc_m; state.ac = ac; state.ac_m = ac_m;
	state.vf = vf; state.af = af;
	state.vm = vm; state.am = am;

	% derivatives
	state.Dx = Dx; state.Dy = Dy; state.DTD = DTD;

	% rot symmetry
	state.Rn = Rn;
end

% DBG, cost function value
if(do_cost)
	cost = cost(1:iter);
end

function y = Ax(x)
xf = x(:, 1:fsize(3)*fsize(4)); xm = x(:, fsize(3)*fsize(4)+1:end);

F(idx_large_f) = xf; M(idx_large_m) = xm;

% forward conv
HF = sum(H.*fft2(F), 4);
if(single_m)
	bHM = b.*real(ifft2(sum(H,4).*fft2(M)));
else
	bHM = b.*real(ifft2(sum(H.*fft2(M), 4)));
end

% transpose conv
yf = real(ifft2(HT.*(HF - fft2(bHM)))); yf = gamma*reshape(yf(idx_large_f),fsize(1)*fsize(2),fsize(3)*fsize(4));
if(single_m)
	ym = real(ifft2(sum(HT,4).*fft2(sum(b.*(bHM - real(ifft2(HF))),3))));
else
	ym = real(ifft2(HT.*fft2(sum(b.*(bHM - real(ifft2(HF))),3)))); 
end
ym = gamma*reshape(ym(idx_large_m), fsize(1)*fsize(2),[]);

% template matching term
%if(any(lambda > 0))
	% multi-template version with formulation lambda*|F_i-T_i|^2 (note: requires const term in the RHS)
	yf = yf + lambda.*xf;

	% original version lambda*|f-m*T|^2: (note: remove const term in the RHS when using this formulation)
	% temp = reshape(T, size(T,1), fsize(3), fsize(4));
	% temp2 = reshape(xm, size(xm,1),1,[]);
	% yf = yf + lambda.*(xf-reshape(temp.*temp2, size(yf,1), []));
	% ym = ym + lambda.*reshape(sum(temp.^2,2).*temp2 - sum(temp.*reshape(xf, size(xf,1), fsize(3), fsize(4)),2), size(ym,1), []);
%end

% cross-img derivative
yf = yf + beta_cross_f*crossDf_DTD(xf);
ym = ym + beta_cross_m*crossDm_DTD(xm);

% rot symmetry
ym = ym + lambda_m0*xm + lambda_R*(Rn*xm);

% regularizers/identity terms
y = [yf ym] + (DTD*x).*beta_tv4 + x*beta_fm;
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

function res = crossD(x, num_channels)
% cross-image forward derivative (img2-img1); expects input format of 'x' as in 'f' and 'm' in the main function
res = x(:,num_channels+1:end)-x(:,1:end-num_channels);
end

function res = crossD_T(x, num_channels)
% transpose of cross-image derivative (img2-img1); expects input format of 'x' as in 'f' and 'm' in the main function but 1 image less (as is output of crossD)
res = [-x(:,1:num_channels) x(:,1:end-num_channels)-x(:,num_channels+1:end) x(:,end-num_channels+1:end)];
end

function res = crossD_DTD(x, num_channels)
% short for crossD_T(crossD(x)) (better memory managenemt when written explicitly)
res = [x(:,1:num_channels)-x(:,num_channels+1:2*num_channels) 2*x(:,num_channels+1:end-num_channels)-x(:,1:end-2*num_channels)-x(:,2*num_channels+1:end) x(:,end-num_channels+1:end)-x(:,end-2*num_channels+1:end-num_channels)];
end
