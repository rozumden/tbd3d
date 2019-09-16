function [x iter relres] = conjgrad(Afun, b, x, tol, maxiter, callback)
% CONJGRAD Direct simple but reasonably fast implementation of conjugate gradient method.
% Major change in 2/2017, use conjgrad_old fro compatibility. Last changed 2/2017.
% Minor change in 8/2017: residual is calculated sloppily as a running residual (may accumulate errors) - correct residual calculation requires additional call to afun, which was removed (this is difference from matlab's pcg)
%
% Afun - function handle: Afun(x) returns A*x (x is not reshaped to column vector before calling Afun, so it can have any shape)
% b - RHS (same shape as 'x')
% x - initial guess, zero if not specified (any shape, result then has the same shape and Afun is called with this shape)
% tol - stopping criterion, L2 norm of the relative residual |Ax-b|/|b| (unsquared) <= tol. Can be [] for fixed number of iterations
% maxiter - stopping criterion, max number of iterations. Can be [] for tolerance only or full-solution.
% callback - function: callback(x_k, iter) called every new update to draw current estimate etc. (iter = iteration number, x_k current estimated in its original shape)
%
% relres - relative residual

if(nargin < 3 || isempty(x))
	x = zeros(size(b));
end
if(nargin < 4 || isempty(tol))
	tol = 0;
end
if(nargin < 5 || isempty(maxiter))
	maxiter = numel(b);
end
if(nargin < 6)
	callback = [];
end

% initialization
b2 = b(:)'*b(:);
tol = tol^2*b2; % tolerance for squared L2 norm of the residual
iter = 0; %r2act = [];
r = b - Afun(x);
p = r;
r2old = r(:)'*r(:); r2new = r2old; % initial residual (btw, this is way faster than sum(r.^2) etc.)

if(r2old > tol)
	% main loop
	for iter=1:maxiter
		% display
		if(~isempty(callback))
			callback(x, iter-1);
		end

		% update
		Ap = Afun(p);
		pAp = p(:)'*Ap(:);
		if(pAp <=0 || ~isnumeric(pAp)) % p can be zero or something - calculate residual and call it a day
			break; 
		end
		alpha = r2old/pAp;
		x = x + alpha*p;

		% residual
		r = r - alpha*Ap; % current residual b-Ax
		r2new = r(:)'*r(:);

		% stopping criterion - residual
		if(r2new <= tol)
			% % actual residual (r2 accumulates roundoff error and can be misleading; requires one more call to afun, not strictly necessary for approximate calculations)
			% temp = Afun(x)-b;
			% r2act = temp(:)'*temp(:);
			% if(r2act <= tol)
			% 	break;
			% end
			% r2act = [];
			break;
		end

		% new direction
		p = r + (r2new/r2old)*p;
		r2old = r2new;
	end
end

if(nargout >= 3)
	% if(isempty(r2act))
	% 	% actual residual (r2 accumulates roundoff error and can be misleading; requires one more call to afun)
	% 	temp = Afun(x)-b;
	% 	r2act = temp(:)'*temp(:);
	% end
	% relres = sqrt(r2act/b2);

	relres = sqrt(r2new/b2);
end

% final display
if(~isempty(callback))
	callback(x, iter);
end
end
