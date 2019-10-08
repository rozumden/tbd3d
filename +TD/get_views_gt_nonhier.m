function [matF, matM, ind] = get_views_gt_nonhier(video, gt_coeffs, m, f, n, params_tbd3d)
if ~exist('params_tbd3d','var')
	params_tbd3d = [];
end
f0_maxiter = IF(isfield(params_tbd3d, 'f0_maxiter'), @()params_tbd3d.f0_maxiter, 30); 
maxiter = IF(isfield(params_tbd3d, 'maxiter'), @()params_tbd3d.maxiter, 20); 


matF = [];
matM = [];
ind = [];
gm = 1;

gap = 0.07;

im = video(:,:,:,1);
video(video < 0) = 0; video(video > 1) = 1;
if size(video,4) < 20
	bgr = median(video(:,:,:, 1:end/2), 4);
else
	bgr = median(video(:,:,:, 1:20), 4);
end
bgr(bgr < 0) = 0; bgr(bgr > 1) = 1;
bgr = bgr.^(gm);
f = f.^(gm);
f = f.*repmat(m,[1 1 3]);

for ci = 1:numel(gt_coeffs)
	coeff = gt_coeffs{ci};
	
	img = video(:,:,:,ci).^(gm);
	
	[h,~,crvlen] = myTrajRender(size2(bgr), coeff, [0 1]);
	h = h / sum(h(:));

	[F0,M0,roi] = estimateFM_motion_template_pw(img, bgr, h, f, m, f, [], 'alpha', 2^-10, 'alpha_m', 2^-12, 'gamma', 1, 'beta_fm', 1e-3, 'lambda', 1e-3, 'maxiter', f0_maxiter, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);

	Ftemplate = repmat(F0, [1 1 1 n]);
	Mtemplate = repmat(M0, [1 1 1 n]);

	Hs = [];
	for k = 1:n
		H = evaluate_vcoeff(size2(im), coeff, [( ( (k-1+gap) )/n) (k/n)]);
		H = double(H);
		H = H ./ sum(H(:));
		Hs = cat(3, Hs, H);
	end
	Hs = Hs / sum(Hs(:));
		
	[Fs,Ms,roi] = estimateFM_motion_template_pw(img, bgr, Hs, Ftemplate, Mtemplate, [], [],'alpha', 3^-10, 'alpha_m', 2^-12,  'gamma', 1, 'beta_fm', 1e-3, 'lambda', 1e-3, 'maxiter', maxiter, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);
	Fs = Fs.^(1/gm); Ms = Ms.^(1/gm);

	matF = cat(4, matF, Fs);
	matM = cat(4, matM, Ms);
	ivs = linspace(ci,ci+1,n+1); 
	ind = [ind ivs(1:end-1)];
end
