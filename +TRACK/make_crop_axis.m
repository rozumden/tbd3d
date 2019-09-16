function [im_c, bgr_c, mbb] = make_crop_axis(im, bgr, ext, pts)
intmethod = 'cubic'; % linear cubic
extvalue = 0;

T = zeros(size(im,1), size(im,2));
pts0 = round(pts);
pts0(pts0 < 1) = 1;
szf = [size(im,2) size(im,1)];
szfm = repmat(szf, size(pts,1), 1);
pts0(pts0 > szf) = szfm(pts0 > szf);
T(sub2ind(size(T), pts0(:,2), pts0(:,1))) = 1;
B = APP.gen_ball(ceil(ext));
T = conv2(T, B, 'same');
[y,x] = find(T>0);
mbb = minBoundingBox([x y]');
% mbbv = mbb(:,[1:end 1]); plot(mbbv(1,:), mbbv(2,:), '-r');
w = norm(mbb(:,1) - mbb(:,2));
h = norm(mbb(:,3) - mbb(:,2));

%% Normalization 1: the image should be blurred in width direction
if w < h 
	mbb = mbb(:,[2 3 4 1]);
	w1 = w;
	w = h;
	h = w1;
end

%% Normalization 2: the image should be from left to right
dist1 = norm(mbb(:,1)-pts(1,:)');
dist2 = norm(mbb(:,2)-pts(1,:)');
if dist2 < dist1
	mbb = mbb(:,[2 1 4 3]);
end

%% Normalization 3: the image should be gravity-friendly 
if mbb(2,4) < mbb(2,1)
	mbb = mbb(:,[4 3 2 1]);
end

ww = ceil(w);
hh = ceil(h);

im_c = zeros(hh,ww,3);
bgr_c = zeros(hh,ww,3);

wdir = mbb(:,2) - mbb(:,1);
wdir = wdir / norm(wdir);
hdir = mbb(:,3) - mbb(:,2);
hdir = hdir / norm(hdir);

wp = mbb(:,1) + [0:(ww-1)].*wdir;
ap = repmat(wp, [1 hh]);
shf = repmat([0:(hh-1)], [ww 1]);
hp = ap + [shf(:)]'.*hdir;

[xx yy] = meshgrid([1:size(im,2)],[1:size(im,1)]);
for k = 1:size(im,3)
	im_c(:,:,k) = (reshape(interp2(xx,yy,double(im(:,:,k)),hp(1,:),hp(2,:),intmethod,extvalue), [ww hh]))';
	bgr_c(:,:,k) = (reshape(interp2(xx,yy,double(bgr(:,:,k)),hp(1,:),hp(2,:),intmethod,extvalue), [ww hh]))';
end