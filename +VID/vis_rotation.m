function [out] = vis_rotation(bgr, Hs, Fs, Ms, gm)
if nargin < 5
	gm = 1/2.2;
end
out = [];
H = Hs;
Fs(Fs < 0) = 0;
for k = 1:size(Fs,4)
	F = Fs(:,:,:,k);
	M = Ms(:,:,:,k);
	if size(Hs, 4) > 1,
		H = Hs(:,:,:,k);
	end
	Mused = M;
	M3 = repmat(Mused, [1 1 3]);
	F(F > M3) = M3(F > M3);

	MH = repmat(conv2(H, Mused, 'same'), [1 1 3]);
	MH(MH > 1) = 1;
	FH = conv2(H, F(:,:,1), 'same');
	FH(:,:,2) = conv2(H, F(:,:,2), 'same');
	FH(:,:,3) = conv2(H, F(:,:,3), 'same');

	im = (bgr.^(1/gm)) .* (1 - MH) + FH; 

	out = cat(4,out,im.^(gm));
end

[cy,cx] = find(H == max(H(:))); 
cx = cx(1); cy = cy(1);

sz = size(M,1);
fc = 5;
% out = out(cy-fc*sz:cy+fc*sz, cx-fc*sz:cx+fc*sz, :, :);

for k = 1:3
	out = cat(4,out,out);
end

% format_file = 'Archival';
format_file = 'Motion JPEG AVI';

out2 = uint8(round(255*out));
if size(Hs, 4) > 1
	VID.write_video(out2,['~/projects/vis/videos/rotation_translation'], format_file, 10);
else
	VID.write_video(out2,['~/projects/vis/videos/rotation'], format_file, 10);
end