function [h f m iter] = blindLoop_FMH_TbD(img, bg, f, m, h, hmask, template, varargin)
% blind alternating loop FM-estimation/H-estimation, either any combination of FM and H can be specified as input

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 100); % max number of f/h alternations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 5e-3); % relative tolerance for h
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % true = debug mode

do_m = IF(isfield(params, 'do_m'), @()params.do_m, true); % false = mask is fixed and not estimated
alpha_f = IF(isfield(params, 'alpha_f'), @()params.alpha_f, 2^-10);
beta_f = IF(isfield(params, 'beta_f'), @()params.beta_f, 2e1*alpha_f);
alpha_ml1 = IF(isfield(params, 'alpha_ml1'), @()params.alpha_ml1, 0);
beta_ml1 = IF(isfield(params, 'beta_ml1'), @()params.beta_ml1, 1e-2*alpha_ml1);
alpha_h = IF(isfield(params, 'alpha_h'), @()params.alpha_h, 2);
beta_h = IF(isfield(params, 'beta_h'), @()params.beta_h, 1e3);
beta_fm = IF(isfield(params, 'beta_fm'), @()params.beta_fm, 1e-3);
lambda = IF(isfield(params, 'lambda'), @()params.lambda, IF(~isempty(template), 1e-1, 0)); % template weight term
lambda_m0 = IF(isfield(params, 'lambda_m0'), @()params.lambda_m0, 0); % mask l2 term |m-m0|

cg_maxiter = IF(isfield(params, 'cg_maxiter'), @()params.cg_maxiter, 25);
cg_tol = IF(isfield(params, 'cg_tol'), @()params.cg_tol, 1e-6);

% hard-coded params
gamma = 1;
fig_h = []; fig_f = [];
if(verbose)
	fig_f = 1; fig_h = 2;
end

% init
fstep_first = isempty(f);
if(isempty(f))
	if(isempty(m))
		error('either `f` or `m` must be specified to determine size');
	end
	f = ones([size(m) size(img,3)]);
elseif(isempty(m))
	m = ones(size2(f));
end
if(isempty(h))
	if(isempty(hmask))
		h = ones(size2(img)); h = h/numel(h);
	else
		h = double(hmask)/nnz(hmask);
	end
end
%if(~do_m) lambda = max(lambda(:)); end % no spatially-variant weighting for fixed mask
m0 = m;
state_h = []; state_f = [];

% main loop
for iter = 1:maxiter
	h_old = h;

	% F-estimation
	if(iter > 1 || fstep_first)
		if(do_m)
			[f m state_f] = estimateFM_motion_template(img, bg, h, f, m, template, state_f, 'alpha', alpha_f, 'beta_tv', beta_f, 'gamma', gamma(min(iter,end)), 'beta_fm', beta_fm, 'alpha_ml1', alpha_ml1, 'beta_ml1', beta_ml1, 'lambda', lambda, 'lambda_m0', lambda_m0, 'm0', m0, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter, 'cg_tol', cg_tol, 'fig', fig_f, 'verbose', verbose, 'outer_iter', iter);
		else
			[f state_f] = estimateF_motion_template(img, bg, h, m, f, template, state_f, 'alpha', alpha_f, 'beta_tv', beta_f, 'gamma', gamma(min(iter,end)), 'beta_f', beta_fm, 'lambda', lambda, 'hard_m', false, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter, 'cg_tol', 1e-6, 'fig', fig_f, 'verbose', verbose, 'outer_iter', iter);
		end
	end
	
	% H-estimation
	[h state_h] = estimateH_motion(img, bg, f, m, h, hmask, state_h, 'alpha', alpha_h, 'beta_lp', beta_h, 'beta_pos', 0, 'gamma', gamma(min(iter,end)), 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter, 'cg_tol', cg_tol, 'fig', fig_h, 'verbose', verbose);
	
	% convergence test
	reldiff2 = sum((h_old(:)-h(:)).^2)/sum(h(:).^2);
	if(verbose)
		fprintf('outer iter=%d, reldiff_h=%.1e\n', iter, sqrt(reldiff2));
	end

	if(reldiff2 < rel_tol^2) % convergence
		break;
	end
end

h(h < 0) = 0;
end