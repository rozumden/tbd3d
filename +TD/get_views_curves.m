function [matF, matM, ind] = get_views_curves(video, curves, m, f)
matF = [];
matM = [];
ind = [];
expos = 0.93;
gm = 1;

gap = 0.07;
mdiff_th = 0.2;
pixels_per_view = 5;
max_n = 8;

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
for ci = 1:numel(curves)
	crv = curves(ci);
	if crv.fit_iv == 0, continue; end
	ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
	nav = max(1, min(max_n, 2^round(log2( round(crv.len / (pixels_per_view * ilen))))));
	nn = nav*ilen;
	n = nav;
	
	for ti = crv.fit_iv(1):crv.fit_iv(2)
		img = video(:,:,:,ti).^(gm);
		if strcmp(crv.type, 'bounce')
			ti_0 = ti - crv.fit_iv(1); crvlen = [];
			h = evaluate_vcoeff(size2(im), crv.coeff, [ti_0/ilen (ti_0+expos)/ilen]);
		else
			n = [];
			[h,~,crvlen] = myTrajRender(size2(bgr), crv.coeff, [ti ti+expos]);
		end
		h = h / sum(h(:));

		[F0,M0,roi] = estimateFM_motion_template_pw(img, bgr, h, f, m, [], [], 'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 30, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
	
		if isempty(n)
			n = max(1, min(max_n, 2^round(log2( round(crvlen / pixels_per_view) ))));
		end

		Fs = F0; Ms = M0;
		HR = cat(4, Fs, repmat(BL,[1 1 1 n-1]));
		for powi = 1:log2(n)
			ni = 2^powi;
			Ftemplate = Fs(:,:,:,repelem(1:size(Fs,4),2));
			Mtemplate = Ms(:,:,:,repelem(1:size(Fs,4),2));

			Hs = [];
			for k = 1:ni
				if strcmp(crv.type, 'bounce')
					kk = k + ni*(ti - crv.fit_iv(1));
					H = evaluate_vcoeff(size2(im), crv.coeff, [( ( (kk-1+gap) )/nn) (kk/nn)]);
				else
					H = myTrajRender(size2(im), crv.coeff, [ti+( (k-1+gap)/ni)*expos ti+((k)/ni)*expos]);
				end
				H = double(H);
				if sum(H(:)) == 0, error('Sum cannot be zero'); end
				H = H ./ sum(H(:));
				Hs = cat(3, Hs, H);
			end
			Hs = Hs / sum(Hs(:));
			
			[Fs,Ms,roi] = estimateFM_motion_template_pw(img, bgr, Hs, Ftemplate, Mtemplate, [], [],'alpha', 3^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-3, 'lambda_m0', 0, 'm0', m, 'maxiter', 20, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);
			Fs = Fs.^(1/gm); Ms = Ms.^(1/gm);
			% Fs = Fs.*m; Ms = Ms.*m;
			
			for ttt = 1:ni, HR = cat(4, HR, Fs(:,:,:,ttt), repmat(BL,[1 1 1 n/ni-1])); end

		end

		matF = cat(4, matF, Fs);
		matM = cat(4, matM, Ms);
		ind = [ind repmat(ti, 1, size(Fs,4))];

		if false
			I1 = montage(HR, 'Size', [size(HR,4)/n n]); C1 = I1.CData;
			I2 = montage(HR0, 'Size', [size(HR0,4)/n n]); C2 = I2.CData;
			I3 = montage(HR00, 'Size', [size(HR00,4)/n n]); C3 = I3.CData;
			HR2 = cat(1, cat(2,f,repmat(BL,[1 3*n-1])), cat(2,C3,C2,C1));
			Fsvis = montage(Fs, 'Size', [1 n]);
			FMsvis = montage(cat(4,repmat(Ms,[1 1 3]),Fs), 'Size', [2 n]);
			im_c = img(roi(1):roi(2), roi(3):roi(4), :);
			bgr_c = bgr(roi(1):roi(2), roi(3):roi(4), :);
			Hs_c = Hs(roi(1):roi(2), roi(3):roi(4), :, :);
			[out] = VID.vis_rotation(bgr_c, Hs_c, Fs, repmat(Mused,[1 1 1 n]), 1);
			[out] = VID.vis_rotation(bgr_c, Hs_c(:,:,:,end/2), Fs, repmat(M,[1 1 1 n]), 1);
		end
	end

end

if false
	inl = logical(ones(1,size(matF,4)));
	for k = 1:size(matF,4)
		F = matF(:,:,:,k);
		M = matM(:,:,:,k);
		if mean(abs(M(m == 1) - m(m == 1))) > mdiff_th
			inl(k) = 0;
		end
	end
	matF = matF(:,:,:,inl);
	matM = matM(:,:,:,inl);
	ind = ind(inl);
end

ntemp = ceil(sqrt(size(matF,4)));
temp = montage(matF, 'Size', [ntemp ntemp]);
reshape([ind zeros(1,ntemp^2-size(matF,4))], [ntemp ntemp])';
