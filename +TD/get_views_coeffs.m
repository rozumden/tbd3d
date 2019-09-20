function [Fs, Ms] = get_views_coeffs(frm, coeff)
if nargin < 2
	[coeff, len, pnts] = postprocc(frm.coeff, [], 100, [0 1]);
end

pixels_per_view = 5;
max_n = 8;

m = double(diskMask([],round(1.5*frm.Radius)));
BL = zeros([size(m) 1]);

n = max(1, min(max_n, 2^round(log2( round(frm.Length / pixels_per_view)))));

[h,~,crvlen] = myTrajRender(size2(frm.im_c), coeff, [0 1]);
h = h / sum(h(:));

[F0,M0,roi] = estimateFM_motion_template_pw(frm.im_c, frm.bgr_c, h, [], m, [], [], 'alpha', 2^-10, 'alpha_m', 3^-9, 'lambda', 1e-1, 'maxiter', 30, 'cg_maxiter', 25, 'cg_tol', 1e-6);

Fs = F0; Ms = M0;
HR = cat(4, Ms, repmat(BL,[1 1 1 n-1]));
for powi = 1:log2(n)
	ni = 2^powi;
	Ftemplate = Fs(:,:,:,repelem(1:size(Fs,4),2));
	Mtemplate = Ms(:,:,:,repelem(1:size(Fs,4),2));

	Hs = [];
	for k = 1:ni
		H = myTrajRender(size2(frm.im_c), coeff, [(k-1)/ni k/ni]);
		H = double(H);
		if sum(H(:)) == 0, error('Sum cannot be zero'); end
		H = H ./ sum(H(:));
		Hs = cat(3, Hs, H);
	end
	Hs = Hs / sum(Hs(:));

	[Fs,Ms,roi] = estimateFM_motion_template_pw(frm.im_c, frm.bgr_c, Hs, Ftemplate, Mtemplate, [], [], ...
		'alpha', 3^-10, 'alpha_m', 3^-9, 'lambda', 1e-3, 'maxiter', 20, 'cg_maxiter', 50, 'cg_tol', 1e-6);
	
	for ttt = 1:ni, HR = cat(4, HR, Ms(:,:,:,ttt), repmat(BL,[1 1 1 n/ni-1])); end

end

% montage(Fs, 'Size', [1 size(Fs,4)]);
montage(Ms, 'Size', [1 size(Fs,4)]);

if false
	[Fs,Ms,roi] = estimateFM_motion_template_pw(frm.im_c, frm.bgr_c, Hs, [], zeros(size(Ms)), [], [], ...
		'alpha', 3^-10, 'alpha_m', 3^-9, 'lambda', 1e-3, 'maxiter', 20, 'cg_maxiter', 50, 'cg_tol', 1e-6);
	
	finit = repmat({ones(size(F0))},1,ni); % initial segments of H; no update of H in convsparseCG_v1z
	minit = repmat({0},1,ni); % independent masks
	[~, Fs, Ms] = convsparseCG_v1z(frm.im_c,frm.bgr_c,finit,minit,Hs);
end
