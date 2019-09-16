function [out] = remove_fmo(video, curves, r)
if size(video,4) < 20
	bgr = median(video(:,:,:, 1:end/2), 4);
else
	bgr = median(video(:,:,:, 1:20), 4);
end
M = double(diskMask([], r*1.3));
out = video;

for ci = 1:numel(curves)
	crv = curves(ci);
	if crv.fit_iv == 0, continue; end
	
	for ti = crv.fit_iv(1):crv.fit_iv(2)
		bgr0 = bgr;
		if ti >= 5
			bgr0 = median(video(:,:,:, (ti-4):ti), 4);
		end
		im = video(:,:,:,ti);

		if strcmp(crv.type, 'bounce')
			BW = myTrajRender(size2(im), crv.coeff, [0 1]);
		else
			BW = myTrajRender(size2(im), crv.coeff, [ti ti+1]);
		end
		BW2 = conv2(BW, M, 'same');
		mask = repmat(BW2,[1 1 3]);
		im(find(mask)) = bgr(find(mask)); 
		out(:,:,:,ti) = im;
	end

end