function [out] = put_fmo(video, curves, M, F, n)
out = video;
FM = F.*repmat(M,[1 1 3]);
for ci = 1:numel(curves)
	crv = curves(ci);
	if crv.fit_iv == 0, continue; end
	
	for ti = crv.fit_iv(1):crv.fit_iv(2)
		for k = 1:n
			ind = n*(ti-1) + k;
			if ind > size(video,4), continue; end
			im = video(:,:,:,ind);

			if strcmp(crv.type, 'bounce')
				ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
				kk = k + n*(ti - crv.fit_iv(1));
				nn = n*ilen;
				H = evaluate_vcoeff(size2(im), crv.coeff, [( (kk-1)/nn) (kk/nn)]);
			else
				H = myTrajRender(size2(im), crv.coeff, [ti+( (k-1)/n) ti+((k)/n)]);
			end
			H = double(H);
			H = H ./ sum(H(:));

			MH = repmat(conv2(H, M, 'same'), [1 1 3]);
			FH = conv2(H, FM(:,:,1), 'same');
			FH(:,:,2) = conv2(H, FM(:,:,2), 'same');
			FH(:,:,3) = conv2(H, FM(:,:,3), 'same');

			im = im .* (1 - MH) + FH; 

			out(:,:,:,ind) = im;
		end
	end

end
		