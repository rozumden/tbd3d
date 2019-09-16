function [T, coeff, s, len] = fitting(h, th, typ, params)
if nargin < 2
	th = 0.4;
end
if nargin < 3
	typ = 1;
end
hext = conv2(h, ones(3), 'same');
mval = max(hext(:));
if mval > th
	[T, coeff, s, len] = FIT.getpoint(h, hext, th, mval);
	return;
end

hm = h/max(h(:));
th1 = 0.1;
th2 = prctile(hm(hm>th1), 50);
hn = h;
hmb = hm.^ (1+th2) > th1;
cnts = conv2(double(hmb), ones(5), 'same');
hn(cnts <= 2) = 0;

if typ == 2
	[T, coeff] = psffit2( (hn/max(hn(:))).^ (1+th2), [th2 th1]); 
elseif typ == 1
	% [T, coeff] = psffit(h, [65 25]/255, 1, 1, 4, 4.3);
	[T, coeff] = psffit(h, params.fit_hard_thresh, params.fit_mass_pcent, ...
			params.fit_num_sampling_cycles, params.fit_inlier_thresh, ...
			params.fit_weighting_sigma, params.fit_max_merge_gap, ...
			params.fit_max_link_gap, params.fit_max_cutoff, ...
			params.fit_max_cos);
else
	[T, coeff] = psffit((hn/max(hn(:))).^ (1+th2), [th2 th1], 1, 1, 4, 4.3);
end
if isempty(coeff)
	[T, coeff, s, len] = DET.getpoint(h, hext, th, mval);
end
[T,~, len] = myTrajRender(size(h), cellfun(@(x)fliplr(x).', coeff, 'UniformOutput', false));

T0s = imdilate(T>0,strel('square',2));
s = sum(h(T0s));
T = T / sum(T(:));

