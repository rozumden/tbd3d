function [out] = put_fmo_pw(video, video0, curves, m, f0, n, do_postprocess)
if ~exist('do_postprocess','var')
	do_postprocess = true;
end

expos = 0.93;
gm = 1;
mdiff_th = 0.2;

im = video(:,:,:,1);
video(video < 0) = 0; video(video > 1) = 1;
video0(video0 < 0) = 0; video0(video0 > 1) = 1;
if size(video0,4) < 20
	bgr = median(video0(:,:,:, 1:end/2), 4);
else
	bgr = median(video0(:,:,:, 1:20), 4);
end
bgr(bgr < 0) = 0; bgr(bgr > 1) = 1;
bgr = bgr.^(gm);
f = f0.^(gm);
f = f.*repmat(m,[1 1 3]);

out = video;

Flast = f; Mlast = m;

for ci = 1:numel(curves)
	crv = curves(ci);
	if crv.fit_iv == 0, continue; end
	ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
	nn = n*ilen;
	
	for ti = crv.fit_iv(1):crv.fit_iv(2)
		img = video0(:,:,:,ti).^(gm);
		if strcmp(crv.type, 'bounce')
			ti_0 = ti - crv.fit_iv(1);
			h = evaluate_vcoeff(size2(im), crv.coeff, [ti_0/ilen (ti_0+1)/ilen]);
		else
			h = myTrajRender(size2(bgr), crv.coeff, [ti ti+expos]);
		end
		h = h / sum(h(:));

		[F0,M0] = estimateFM_motion_template(img, bgr, h, f, m, f, [], 'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
		if do_postprocess
			[F0,M0] = postprocessFM(F0,M0);
			% F0 = F0.*M0;
		end
	
		Fdiscrete = F0; Mdiscrete = M0;
		if mean(abs(M0(m == 1) - m(m == 1))) > mdiff_th
			Fdiscrete = f; Mdiscrete = m;
		end

		Hs = [];
		for k = 1:n
			ind = n*(ti-1) + k;
			if ind > size(video,4), continue; end
			im = video(:,:,:,ind);

			if strcmp(crv.type, 'bounce')
				kk = k + n*(ti - crv.fit_iv(1));
				H = evaluate_vcoeff(size2(im), crv.coeff, [( (kk-1)/nn) (kk/nn)]);
			else
				H = myTrajRender(size2(im), crv.coeff, [ti+( (k-1)/n) ti+((k)/n)]);
			end
			H = double(H);
			if sum(H(:)) == 0, continue; end
			H = H ./ sum(H(:));
			% H = H / n;
			Hs = cat(3, Hs, H);
		end
		if isempty(Hs)
			continue
		end

		Hs = Hs / sum(Hs(:));
		% [Fs,Ms,roi] = estimateFM_motion_pw(img, bgr, Hs, repmat(Fdiscrete,[1 1 1 n]), Mdiscrete, [],'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
		[Fs,Ms,roi] = estimateFM_motion_pw(img, bgr, Hs, repmat(Fdiscrete,[1 1 1 n]), Mdiscrete ,  [],'alpha', 3^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 10, 'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6);
		
		Fs0 = []; Hs0 = [];
		for k = 1:n
			ind = n*(ti-1) + k;
			if ind > size(video,4), continue; end
			im = video(:,:,:,ind);

			F = Fs(:,:,:,k);
			M = Ms;
			H = Hs(:,:,k);
			if do_postprocess
				[F,M] = postprocessFM(F,M);
				F = F.*M;
			end

			F = F.^(1/gm);
			M = M.^(1/gm);

			if mean(abs(M(:) - m(:))) > mdiff_th && ~isempty(Flast)
				F = Flast; M = Mlast;
			end
			Flast = F; Mlast = M;

			Mused = M;
			Mused(m == 0) = 0;
			M3 = repmat(Mused, [1 1 3]);
			F(F > M3) = M3(F > M3);

			H = H / sum(H(:));
			MH = repmat(conv2(H, Mused, 'same'), [1 1 3]);
			MH(MH > 1) = 1;
			FH = conv2(H, F(:,:,1), 'same');
			FH(:,:,2) = conv2(H, F(:,:,2), 'same');
			FH(:,:,3) = conv2(H, F(:,:,3), 'same');

			im = im .* (1 - MH) + FH; 

			% if k == 3 && ti > 15, keyboard; end

			out(:,:,:,ind) = im;

			Fs0 = cat(4, Fs0, F); Hs0 = cat(4,Hs0,H);
		end
		if false
			% Fsvis = montage(Fs0, 'Size', [1 n]);
			% FMsvis = montage(cat(4,repmat(Ms0,[1 1 3]),Fs0), 'Size', [2 n]);
			im_c = img(roi(1):roi(2), roi(3):roi(4), :);
			bgr_c = bgr(roi(1):roi(2), roi(3):roi(4), :);
			Hs_c = Hs0(roi(1):roi(2), roi(3):roi(4), :, :);
			[out] = VID.vis_rotation(bgr_c, Hs_c, Fs, repmat(Mused,[1 1 1 n]), 1);
			[out] = VID.vis_rotation(bgr_c, Hs_c(:,:,:,end/2), Fs, repmat(M,[1 1 1 n]), 1);
		end
		if ti > 18, keyboard; end
	end

end

		