function [matF, matM, ind] = get_views_hs(video, curves, par, M, nf, use_gt)
if ~exist('use_gt','var')
	use_gt = true;
end

matF = [];
matM = [];
ind = [];
expos = 0.93;
gap = 0.07;
dif_th = 0.1;

bgr = median(video(:,:,:, 1:20), 4);
im = video(:,:,:,1);
M3 = repmat(M,[1 1 3]) > 0;
if use_gt
	for k = 1:numel(par)
		for kk = 1:size(par{k},2)
			ind = nf*(k-1) + kk;
			img = video(:,:,:,ind);
			pr = par{k}(:,kk);

			F = zeros([size(M) 3]);
			H = zeros(size2(img));
			H(pr(2),pr(1),:) = 1;
			H = conv2(H,M,'same') > 0;
			H3 = repmat(H,[1 1 3]);
			F(M3) = img(H3);
			matF = cat(4, matF, F);
		end
	end
else 
	for ci = 1:numel(curves)
		crv = curves(ci);
		if crv.fit_iv == 0, continue; end
		% if strcmp(crv.type,'prediction'), continue; end
		ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
		nn = nf*ilen;
		
		for ti = crv.fit_iv(1):crv.fit_iv(2)
			for k = 1:nf
				img = video(:,:,:, nf*(ti-1) + k);
				if strcmp(crv.type, 'bounce')
					kk = k + nf*(ti - crv.fit_iv(1));
					H = evaluate_vcoeff(size2(im), crv.coeff, [( ( (kk-1+gap) )/nn) (kk/nn)]);
				else
					H = myTrajRender(size2(im), crv.coeff, [ti+( (k-1+gap)/nf)*expos ti+((k)/nf)*expos]);
				end
				H = double(H);
				H = H ./ sum(H(:));
				[y,x] = find(H>0);
				x = round(mean(x)); y = round(mean(y));

				F = zeros([size(M) 3]);
				H = zeros(size2(img));
				H(y,x,:) = 1;
				H = conv2(H,M,'same') > 0;
				H3 = repmat(H,[1 1 3]);
				F(M3) = img(H3);

				matF = cat(4, matF, F);
			end
		end
	end
end

