function [f_img state cost_data cost_tv] = estimateF_motion_template(g, b, h, mask, f, T, state, varargin)
% FIXME
%
% g,b - input, background resp., double RGB, same size
% h - psf, same 2-size as g,b
% mask - logical mask, same 2-size as f
% state - strcut containing auxiliary varibles required to resume optimization; simply 'state' as in output

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 50); % max number of outer iterations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 2e-3); % relative between iterations difference for outer admm loop
cg_maxiter = IF(isfield(params, 'cg_maxiter'), @()params.cg_maxiter, 50); % max number of inner CG iterations
cg_tol = IF(isfield(params, 'cg_tol'), @()params.cg_tol, 1e-5); % tolerance for relative residual of inner CG iterations ('h' subproblem); can be several values which are used sequentially each time the convergence criterion holds
gamma = IF(isfield(params, 'gamma'), @()params.gamma, 1); % data term weight
alpha = IF(isfield(params, 'alpha'), @()params.alpha, 2^-10); % TV regularizer weight (f)
lambda = IF(isfield(params, 'lambda'), @()params.lambda, 0); % template L2 term weight
beta_tv = IF(isfield(params, 'beta_tv'), @()params.beta_tv, 10*alpha); % splitting vx/vy=Df due to TV regularizer
beta_f = IF(isfield(params, 'beta_f'), @()params.beta_f, 1e-3); % splitting vf=f due to masking+positivity constraint, penalty weight
hard_m = IF(isfield(params, 'hard_m'), @()params.hard_m, false); % if true, masking is forced and 'f' is calculated by cg iterations; if 0 masking is enforced by add penalty and f subproblem is solved without cg (direct inverse in ft)
fig = IF(isfield(params, 'fig'), @()params.fig, []); % figure handle or id for plotting intermediate results; set to [] to disable plotting
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % 1=display info message in each iteration, 0=don't
do_cost = IF(isfield(params, 'do_cost'), @()params.do_cost, false); % true = calculate cost function value (and related quantities) in each iteration (takes time ~ 1CG iteration); the result it returned in 'cost' struct; false = don't
outer_iter = IF(isfield(params, 'outer_iter'), @()params.outer_iter, 1); % iter in the outer loop, for dbg purposes

if(verbose)
	timer = tic(); % DBG
end

rel_tol = rel_tol^2;

% crop images to minimal relevant roi (speeds-up the calculation of FFTs)
roi = boundingBox(h);
fsize = [size(mask) size(g,3)];
radius = ceil((fsize(1:2)-1)/2);
roi = roi + [-radius(1), radius(1), -radius(2), radius(2)]; % relevant region in the input image
% pad_lt = ceil((fsize(1:2)-1)/2); pad_br = floor((fsize(1:2)-1)/2);
% roi = roi + [-pad_lt(1), +pad_br(1), -pad_lt(2), +pad_br(2)]; % relevant region in the input image
roi = min(max(roi, [1 1 1 1]), size2(h,[1 1 2 2]));
g = g(roi(1):roi(2), roi(3):roi(4), :);
b = b(roi(1):roi(2), roi(3):roi(4), :);
h = h(roi(1):roi(2), roi(3):roi(4));
gsize = size2(g, 1:3);

% unknowns
mask3 = repmat(mask, [1 1 fsize(3)]);
if(hard_m)
	if(isempty(f))
		f = zeros(nnz(mask), fsize(3)); % RGB channels of 'f' a column vectors
	else
		f = vec3(f); % RGB channels of 'f' a column vectors
		f = f(mask, :);
	end
	idx_large = psfshift_idx(mask3, gsize); % indices for quick conversion between small 'f' (actual 'f') and 'F' (psf-shifted large version of 'f' for conv in FT)
	if(lambda > 0)
		T = vec3(T); T = T(mask, :); % template, is any
	else
		T = 0;
	end

	% derivatives
	[Dx Dy] = createDerivatives(mask);
	DTD = Dx.'*Dx + Dy.'*Dy;
else % calculation in FT
	if(isempty(f))
		f = zeros(prod(fsize),1); % all channels together
	else
		f = f(:); % all channels together
	end
	idx_large = psfshift_idx(fsize, gsize); % indices for quick conversion between small 'f' (actual 'f') and 'F' (psf-shifted large version of 'f' for conv in FT)
	if(lambda > 0)
		FT = zeros(gsize); FT(idx_large) = T.*mask; % template
	else
		FT = 0;
	end

	% preallocation for FT calculations
	VX = zeros(gsize); AX = zeros(gsize);
	VY = zeros(gsize); AY = zeros(gsize);
	VF = zeros(gsize); AF = zeros(gsize);
	
	% derivatives in FT
	Dx = fft2(psfshift([1 -1], gsize(1:2))); Dy = fft2(psfshift([1; -1], gsize(1:2)));
	DTD = beta_tv*(conj(Dx).*Dx + conj(Dy).*Dy);

	cg_iter = 0; cg_residual = 0;
end

F = zeros(gsize); % prealocation for FT calculation (large psfshifted versions of f/m for conv in FT)

% initial vals
if(nargin >= 6 && ~isempty(state))
	% retrive stored variables
	vx = state.vx; vy = state.vy; ax = state.ax; ay = state.ay; % v=Df splitting due to TV and its assoc. Lagr. mult.
	vf = state.vf; af = state.af; % vf=f splitting due to positivity and f=0 outside mask constraint
else
	vx = zeros(size(f)); vy = zeros(size(f)); ax = 0; ay = 0; % v=Df splitting due to TV and its assoc. Lagr. mult.
	vf = 0; af = 0; % vf=f splitting due to positivity and f=0 outside mask constraint
end

% precompute FT
H = fft2(h); HT = conj(H); HTH = gamma*HT.*H;

% precompute const RHS/LHS for 'f' subproblem
rhs_0 = g-b + b.*real(ifft2(H.*fft2(psfshift(mask, gsize(1:2)))));
if(hard_m)
	rhs_const = real(ifft2(HT.*fft2(rhs_0)));
	rhs_const = gamma*reshape(rhs_const(idx_large), [], fsize(3)) + lambda*T;
else
	rhs_const = gamma*HT.*fft2(rhs_0) + lambda*fft2(FT); % leave in the FT
end

% DBG - cost function evolution
cost_data = zeros(1,maxiter);
cost_tv = zeros(1,maxiter);

% derivatives (note: remove in release version and move the other derivative calculation to the beginning of admm loop)
if(hard_m)
	fdx = Dx*f; fdy = Dy*f;
else
	F(idx_large) = f;
	FF = fft2(F);
	fdx = real(ifft2(Dx.*FF)); fdx = fdx(idx_large);
	fdy = real(ifft2(Dy.*FF)); fdy = fdy(idx_large);
end

% ADMM loop
for iter = 1:maxiter
	f_old = f;
	
	% dual/auxiliary var updates updates
	% vx/vy minimization (splitting due to TV regularizer)
	val_x = fdx+ax; val_y = fdy+ay;
	[shrink_factor] = soft_thresh(sqrt(val_x.^2 + val_y.^2), alpha/beta_tv); % isotropic TV
	%[shrink_factor] = soft_thresh(sqrt(sum(val_x.^2,3) + sum(val_y.^2,3)), alpha/beta_tv); % isotropic color TV
	% NOTE: anisotropix TV is factor_x = soft_thresh(abs(val_x), alpha/beta_tv), the same for _y (ie different shrink factors for val_x and val_y)
	vx = val_x.*shrink_factor;
	vy = val_y.*shrink_factor;

	% 'a' step
	ax = ax + fdx - vx;
	ay = ay + fdy - vy;

	% vf minimization (positivity)
	if(beta_f > 0)
		vf = f + af;
		vf(vf < 0) = 0; % infinite penalty for negative 'f' values
		vf(vf > 1) = 1;

		% masking
		if(~hard_m)
			vf(~mask3) = 0;
		end
		
		% 'af' step
		af = af + f - vf;
	end
	
	% f step
	if(hard_m)
		rhs = rhs_const + beta_tv*(Dx.'*(vx-ax)+Dy.'*(vy-ay)) + beta_f*(vf-af);
		[f, cg_iter, cg_residual] = conjgrad(@Ax, rhs, f, cg_tol, cg_maxiter);
	else
		VX(idx_large) = vx; AX(idx_large) = ax;	VY(idx_large) = vy; AY(idx_large) = ay;	VF(idx_large) = vf; AF(idx_large) = af;
		FF = (rhs_const + beta_tv*(conj(Dx).*fft2(VX-AX) + conj(Dy).*fft2(VY-AY)) + beta_f*fft2(VF-AF))./(HTH + DTD + (beta_f + lambda));
		f = real(ifft2(FF)); f = f(idx_large);
	end

	% NOTE: this is where the convergence tests should be in the release version

	% derivatives
	if(hard_m)
		fdx = Dx*f; fdy = Dy*f;
	else
		fdx = real(ifft2(Dx.*FF)); fdx = fdx(idx_large);
		fdy = real(ifft2(Dy.*FF)); fdy = fdy(idx_large);
	end

	% cost (DBG)
	% FIXME: remove in the release version, calculating data term error costs several uneccessary ffts similar to one CG iteration
	if(do_cost)
		if(hard_m)
			F(idx_large) = f; % large psf-shifted version for FT
			HF = real(ifft2(H.*fft2(F)));
		else
			HF = real(ifft2(H.*FF));
		end
		cost_data(iter) = sum(abs(reshape(HF-rhs_0, [], 1)).^2)/2;
		cost_tv(iter) = alpha*sum(sqrt(fdx(:).^2+fdy(:).^2));
		% FIXME: neni cost za template term
	end
	
	% FIXME: move to convergence test section
	df = f-f_old;
	rel_diff2 = (df(:)'*df(:))/(f(:)'*f(:)); % relative l2 norm squared (between-iterations difference)

	% DBG
	if(verbose)
		fprintf('F: iter=%d, cg_iter=%d, cgres=%.1e, reldiff=%.1e, time=%.2fs\n', iter+outer_iter-1, cg_iter, cg_residual, sqrt(rel_diff2), toc(timer));
		timer = tic();
	end

	% plot (DBG)
	if(~isempty(fig))
		figure(fig);
		if(hard_m)
			f_img = zeros(fsize); f_img(mask3) = f;
		else
			f_img = reshape(f, fsize);
		end
		imshow(f_img); title(sprintf('iter=%d, reldiff=%.1e', iter+outer_iter-1, sqrt(rel_diff2)));
		% subplot(223);
		% yyaxis('left'); plot(1:iter, cost_data(1:iter)+cost_tv(1:iter));
		% yyaxis('right'); plot(1:iter, cost_data(1:iter));
		% title(sprintf('cost=%.2e, data=%.2e', cost_data(iter)+cost_tv(iter), cost_data(iter)));
		drawnow;
	end

	% convergence testing
	if(rel_diff2 < rel_tol)
		break;
	end
end

% reshape images
if(hard_m)
	f_img = zeros(fsize); f_img(mask3) = f;
else
	f_img = reshape(f, fsize);
end

% return state
if(nargout >= 2)
	% return auxiliaries
	state = struct;
	state.vx = vx; state.vy = vy; state.ax = ax; state.ay = ay;
	state.vf = vf; state.af = af;
end

% DBG, cost function value
cost_data = cost_data(1:iter);
cost_tv = cost_tv(1:iter);

function y = Ax(x)
F(idx_large) = x; % large psf-shifted version for FT
y = real(ifft2(HTH.*fft2(F)));
y = reshape(y(idx_large), [], fsize(3)) + beta_tv*DTD*x + (beta_f + lambda)*x;
end
end

function [shrink_factor] = soft_thresh(val_norm, amount) % helper function for vector version of soft thresholding for l1 minimization
nz = val_norm > amount;
shrink_factor = zeros(size(val_norm));
shrink_factor(nz) = (val_norm(nz)-amount)./val_norm(nz);
end