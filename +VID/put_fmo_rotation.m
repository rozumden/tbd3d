function [out] = put_fmo_rotation(video, video0, curves, m, f0, n)
subframe_F = true;

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
% f = meanrotate(f);
f = f.*repmat(m,[1 1 3]);

out = video;

Flast = []; Mlast = [];

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

		try
			[F0,M0,roi] = estimateFM_motion_template(img, bgr, h, f, m, f, [], 'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
			[F0,M0] = postprocessFM(F0,M0);
			% F0 = F0.*M0;
		catch
			disp('Error in estimateFM_motion_template');
			F0 = f; 
			M0 = m;
		end
		Fdiscrete = F0; Mdiscrete = M0;
		if mean(abs(M0(m == 1) - m(m == 1))) > mdiff_th
			Fdiscrete = f; Mdiscrete = m;
		end

		Fs = []; Ms = []; Hs = [];
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

			if subframe_F
				Hmult = H / n;
				try
					% [F,M] = estimateFM_motion_template(img, bgr, Hmult, f, m, f, [], 'alpha', 2^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, 'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', 1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 1, 'rel_tol', 0, 'cg_maxiter', 25, 'cg_tol', 1e-6);
					[F,M] = estimateFM_motion_template(img, bgr, Hmult, ...
						 Fdiscrete, Mdiscrete, Fdiscrete, [], ...
						'alpha', 3^-10, 'beta_tv', 2e1*2^-10, 'gamma', 1, ...
						'beta_fm', 1e-3, 'alpha_ml1', 0, 'beta_ml1', 0, 'lambda', ...
						1e-1, 'lambda_m0', 0, 'm0', m, 'maxiter', 10, 'rel_tol', 0, ...
						'cg_maxiter', 50, 'cg_tol', 1e-6);
					[F,M] = postprocessFM(F,M);
					% F = F.*M;
					F = F.^(1/gm);
					M = M.^(1/gm);
					% if k == 3 && ti > 15, keyboard; end
				catch
					disp('Error in estimateFM_motion_template');
					if ~isempty(Flast)
						F = Flast; M = Mlast;
					else
						F = f; M = m;
					end
				end
			end
			if mean(abs(M(:) - m(:))) > mdiff_th && ~isempty(Flast)
				F = Flast; M = Mlast;
			end
			Flast = F; Mlast = M;
			
			Mused = M;
			Mused(m == 0) = 0;
			M3 = repmat(Mused, [1 1 3]);
			F(F > M3) = M3(F > M3);

			MH = repmat(conv2(H, Mused, 'same'), [1 1 3]);
			MH(MH > 1) = 1;
			FH = conv2(H, F(:,:,1), 'same');
			FH(:,:,2) = conv2(H, F(:,:,2), 'same');
			FH(:,:,3) = conv2(H, F(:,:,3), 'same');

			im = im .* (1 - MH) + FH; 

			out(:,:,:,ind) = im;

			Fs = cat(4, Fs, F); Ms = cat(4, Ms, Mused);
			Hs = cat(4, Hs, H);
		end
		if false
			% Fsvis = montage(Fs, 'Size', [1 n]);
			% FMsvis = montage(cat(4,repmat(Ms,[1 1 3]),Fs), 'Size', [2 n]);
			bgr_c = bgr(roi(1):roi(2), roi(3):roi(4), :);
			Hs_c = Hs(roi(1):roi(2), roi(3):roi(4), :, :);
			[out] = VID.vis_rotation(bgr_c, Hs_c, Fs, Ms, 1);
			[out] = VID.vis_rotation(bgr_c, Hs_c(:,:,:,end/2), Fs, Ms, 1);
		end
		% if ti > 18, keyboard; end
	end

end

		