function [Fs,Ms,Hs] = put_fmo_frame(fr, n, gm)
if nargin < 2
	n = 16;
end
if nargin < 3
	gm = 1/2.2;
end

st_t = 0;
en_t = 1;
len_t = en_t - st_t;

img = fr.im_c.^(1/gm);
bgr = fr.bgr_c.^(1/gm);
coeff = fr.coeff;
if isempty(coeff)
	% [params, cfg] = EVAL.get_params(false);
	% hmask = conv2(fr.T, double(diskMask(30)), 'same') > 0;
	% [h, f, m, T, coeff, len, s] = TRACK.TbD_loop(img, bgr, fr.M, [], [], fr.T, hmask, [], [], params, params.s_th, 10, false);
	coeff = fr.coeffi;
end

for c = 1:numel(coeff)
	coeff{c} = fliplr(coeff{c}).';
end

h = evaluate_vcoeff(size2(img), coeff, [st_t en_t]);
h = h / sum(h(:));
M = diskMask(double(round(2*fr.Radius*0.9)));
% M = fr.M;
[f,m] = estimateFM_motion_template(img, bgr, h, [], M, [], [], ...
	'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, ...
	'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 0, 'lambda_m0', 0, ...
	'm0', M, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
% [f,m] = postprocessFM(f,m);
% f = fr.f; 
% m = fr.M;

Fs = [];
Ms = [];
Hs = [];
for k = 1:n
	H = evaluate_vcoeff(size2(img), coeff, [st_t+((k-1)/n)*len_t st_t+(k/n)*len_t]);
	[y,x] = find(H); plot(x,y,'.r'); drawnow;
	H = double(H);
	if sum(H(:)) == 0, continue; end
	H = H ./ sum(H(:));
	Hmult = H / n;
	[F,M] = estimateFM_motion_template(img, bgr, Hmult, f, m, f, [], ...
		'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, ...
		 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, ...
		 'm0', m, 'maxiter', 10, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);
	[F,M] = postprocessFM(F,M);
	Fs = cat(4, Fs, F);
	Ms = cat(4, Ms, M);
	Hs = cat(4, Hs, H);
end

% Fsvis = montage(Fs, 'Size', [5 n/5]);
FMsvis = montage(cat(4,repmat(Ms,[1 1 3]),Fs), 'Size', [2 n]);

