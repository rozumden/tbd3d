function [mat] = img2mat(fldr, pth, st, en)
if ~exist('st','var')
	st = 1;
end
vid = VID.VideoImg(fldr, pth, '');
if ~exist('en','var')
	en = vid.sz;
end
sz = en - st + 1;
mat = zeros(vid.h,vid.w,3,sz);
for k = st:en
	mat(:,:,:,k-st+1) = vid.get_frame(k);
end
