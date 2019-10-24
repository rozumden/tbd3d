function [matF, matM, ind] = get_views_gt(video, gt_coeffs, m, f, n, params_tbd3d)
if ~exist('params_tbd3d','var')
	params_tbd3d = [];
end
f0_maxiter = IF(isfield(params_tbd3d, 'f0_maxiter'), @()params_tbd3d.f0_maxiter, 30); 
maxiter = IF(isfield(params_tbd3d, 'maxiter'), @()params_tbd3d.maxiter, 20); 
lambda_R = IF(isfield(params_tbd3d, 'lambda_R'), @()params_tbd3d.lambda_R, 1e-2); 
alpha_cross_f = IF(isfield(params_tbd3d, 'alpha_cross_f'), @()params_tbd3d.alpha_cross_f, 2^-12); 
alpha_cross_m = IF(isfield(params_tbd3d, 'alpha_cross_m'), @()params_tbd3d.alpha_cross_m, 2^-12); 
do_intervals = IF(isfield(params_tbd3d, 'do_intervals'), @()params_tbd3d.do_intervals, false); 


matF = [];
matM = [];
ind = [];
expos = 0.93;
gm = 1;
max_n = n;

gap = 0.07;
mdiff_th = 0.2;
pixels_per_view = 5;

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
BL = zeros(size(f));
BLM = zeros(size(m));

for ci = 1:numel(gt_coeffs)
	coeff = gt_coeffs{ci};
	
	img = video(:,:,:,ci).^(gm);
	
	[h,~,crvlen] = myTrajRender(size2(bgr), coeff, [0 1]);
	if sum(h(:)) == 0
		continue;
	end
	
	h = h / sum(h(:));

	[F0,M0,roi] = estimateFM_motion_pw(img, bgr, h, f, m, f, [], 'lambda_R', lambda_R, 'alpha_f', 2^-10, 'alpha_m', 2^-12, 'lambda', 1e-3, 'maxiter', f0_maxiter, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);
	% F0 = F0.*m; M0 = M0.*m;

	Fs = F0; Ms = M0;

	if do_intervals
		n = max(1, min(max_n, 2^round(log2( round(crvlen / pixels_per_view)))));
	end

	HR = cat(4, Fs, repmat(BL,[1 1 1 n-1]));
	HRM = cat(3, Ms, repmat(BLM,[1 1 n-1]));
	for powi = 1:log2(n)
		ni = 2^powi;
		Ftemplate = Fs(:,:,:,repelem(1:size(Fs,4),2));
		Mtemplate = Ms(:,:,repelem(1:size(Fs,4),2));

		Hs = [];
		for k = 1:ni
			H = evaluate_vcoeff(size2(im), coeff, [( ( (k-1+gap) )/ni) (k/ni)]);
			H = double(H);
			H = H ./ sum(H(:));
			Hs = cat(3, Hs, H);
		end
		Hs = Hs / sum(Hs(:));
		
		% Ftemplate = repmat(F0,[1 1 1 ni]); Mtemplate = repmat(M0,[1 1 1 ni]);
		% lambda_w = (ni/n/5) * 1e-2; alpha_w = (2+ni/n)^-10;
		lambda_w = 1e-3; alpha_w = 2^-12;
		[Fs,Ms,roi] = estimateFM_motion_pw(img, bgr, Hs, Ftemplate, Mtemplate, [], [],'alpha_f', alpha_w, 'alpha_m', 2^-12,  'lambda', lambda_w, 'maxiter', maxiter, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6, ...
			'alpha_cross_f', alpha_cross_f, 'alpha_cross_m', alpha_cross_m, 'lambda_R', lambda_R);
		Fs = Fs.^(1/gm); Ms = Ms.^(1/gm);
		
		for ttt = 1:ni, HR = cat(4, HR, Fs(:,:,:,ttt), repmat(BL,[1 1 1 n/ni-1])); end
		for ttt = 1:ni, HRM = cat(3, HRM, Ms(:,:,ttt), repmat(BLM,[1 1 n/ni-1])); end
	end
	HRM = permute(repmat(HRM,[1 1 1 1]),[1 2 4 3]);

	% keyboard

	matF = cat(4, matF, Fs);
	matM = cat(3, matM, Ms);
	ivs = linspace(ci,ci+1,n+1); 
	ind = [ind ivs(1:end-1)];
end
