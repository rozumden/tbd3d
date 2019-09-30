function [f_img,m_img] = estimateFM_mc(video0, curves, m, f)
expos = 0.93;
gm = 1;
mdiff_th = 0.2;

video0(video0 < 0) = 0; video0(video0 > 1) = 1;
if size(video0,4) < 20
	bgr = median(video0(:,:,:, 1:end/2), 4);
else
	bgr = median(video0(:,:,:, 1:20), 4);
end
bgr(bgr < 0) = 0; bgr(bgr > 1) = 1;
bgr = bgr.^(gm);
f = f.^(gm);
f = meanrotate(f);
f = f.*repmat(m,[1 1 3]);

g = {};
b = {};
hs = {};
for ci = 1:numel(curves)
	crv = curves(ci);
	if crv.fit_iv == 0, continue; end
	ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
	
	for ti = crv.fit_iv(1):crv.fit_iv(2)
		img = video0(:,:,:,ti).^(gm);
		if strcmp(crv.type, 'bounce')
			ti_0 = ti - crv.fit_iv(1);
			h = evaluate_vcoeff(size2(bgr), crv.coeff, [ti_0/ilen (ti_0+1)/ilen]);
		else
			h = myTrajRender(size2(bgr), crv.coeff, [ti ti+expos]);
		end
		h = h / sum(h(:));
		g = [g {img}];
		b = [b {bgr}];
		hs = [hs {h}];
	end
end

[f_img,m_img] = estimateFM_motion_template_mc(g, b, hs, f, m, f, [],'alpha', 2^-8, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'maxiter', 50, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);

% [F0,M0] = estimateFM_motion_template(img, bgr, h, f, m, f, [], 'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
% [F0,M0] = postprocessFM(F0,M0);
