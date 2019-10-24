function [h f m iter res_per_iter] = blindLoop_FMH_TbD(img, bg, f, m, h, hmask, varargin)
% blind alternating loop FM-estimation/H-estimation, any combination of FM and H can be specified as input

% params
params = IF(length(varargin) == 1 && isstruct(varargin{1}), @()varargin{1}, struct(varargin{:})); % allows passing optional args as a struct or as key-value pairs
maxiter = IF(isfield(params, 'maxiter'), @()params.maxiter, 100); % max number of f/h alternations
rel_tol = IF(isfield(params, 'rel_tol'), @()params.rel_tol, 5e-3); % relative tolerance for h
verbose = IF(isfield(params, 'verbose'), @()params.verbose, false); % true = debug mode
do_m = IF(isfield(params, 'do_m'), @()params.do_m, true); % false = mask is fixed and not estimated
gamma = IF(isfield(params, 'gamma'), @()params.gamma, 1); % data term weight for both FM and H
alpha_f = IF(isfield(params, 'alpha_f'), @()params.alpha_f, 2^-12);
beta_f = IF(isfield(params, 'beta_f'), @()params.beta_f, 10*alpha_f);
alpha_m = IF(isfield(params, 'alpha_m'), @()params.alpha_m, 2^-12); % mask TV weight
beta_m = IF(isfield(params, 'beta_m'), @()params.beta_m, 10*alpha_m);
beta_fm = IF(isfield(params, 'beta_fm'), @()params.beta_fm, 1e-3);
alpha_h = IF(isfield(params, 'alpha_h'), @()params.alpha_h, 0);
beta_h = IF(isfield(params, 'beta_h'), @()params.beta_h, 1e3);
sum1 = IF(isfield(params, 'sum1'), @()params.sum1, true);
lambda = IF(isfield(params, 'lambda'), @()params.lambda, 0); % template weight term
lambda_m0 = IF(isfield(params, 'lambda_m0'), @()params.lambda_m0, 0); % mask l2 term |m-m0|
lambda_R = IF(isfield(params, 'lambda_R'), @()params.lambda_R, 0); % mask rotation symmetry weight term

% hard-coded params
fig_h = []; fig_f = [];
if(verbose)
	fig_f = 1; fig_h = 2;
end
cg_maxiter_fm = 25;
cg_maxiter_h = 25;

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
m0 = m; template = f;
state_h = []; state_f = [];

if(nargout >= 5)
	res_per_iter = struct([]);
end

% main loop
for iter = 1:maxiter
	h_old = h;

	% F-estimation
	if(iter > 1 || fstep_first)
		if(do_m)
			[f m state_f] = estimateFM_motion_template(img, bg, h, f, m, template, state_f, 'alpha_f', alpha_f, 'beta_f', beta_f, 'alpha_m', alpha_m, 'beta_m', beta_m, 'gamma', gamma(min(iter,end)), 'beta_fm', beta_fm, 'lambda', lambda, 'lambda_m0', lambda_m0, 'm0', m0, 'lambda_R', lambda_R, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter_fm, 'cg_tol', 1e-5, 'fig', fig_f, 'verbose', verbose, 'outer_iter', iter);
		else
			[f state_f] = estimateF_motion_template(img, bg, h, m, f, template, state_f, 'alpha_f', alpha_f, 'beta_f', beta_f, 'gamma', gamma(min(iter,end)), 'beta_f', beta_fm, 'lambda', lambda, 'hard_m', false, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter_fm, 'cg_tol', 1e-5, 'fig', fig_f, 'verbose', verbose, 'outer_iter', iter);
		end
	end
	
	% H-estimation
	[h state_h] = estimateH_motion(img, bg, f, m, h, hmask, state_h, 'alpha_h', alpha_h, 'beta_h', beta_h, 'sum1', sum1, 'gamma', gamma(min(iter,end)), 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', cg_maxiter_h, 'cg_tol', 1e-5, 'fig', fig_h, 'verbose', verbose, 'outer_iter', iter);
	
	% convergence test
	reldiff2 = sum((h_old(:)-h(:)).^2)/sum(h(:).^2);
	if(verbose)
		fprintf('outer iter=%d, reldiff_h=%.1e\n', iter, sqrt(reldiff2));
	end

	if(nargout <= 5)
		res_per_iter(iter).f = f;
		res_per_iter(iter).m = m;
		res_per_iter(iter).h = h;
		res_per_iter(iter).reldiff_h = sqrt(reldiff2);
	end

	if(reldiff2 < rel_tol^2) % convergence
		break;
	end
end

h(h < 0) = 0;
end
