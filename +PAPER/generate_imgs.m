function [] = generate_imgs(matF, frame, matF_hs, iv)
if nargin < 4
	iv = 41:0.5:43;
end
n = size(matF,4) / numel(frame);
frmt = '.png';

for i = iv
	F3D = matF(:,:,:,(i-1)*n + 1); 
	GT = matF_hs(:,:,:,(i-1)*n + 1); 
	imwrite(GT,['~/gt_' num2str(i) frmt]);
	imwrite(F3D,['~/tbd3d_' num2str(i) frmt]);
	if round(i) == i
		imwrite(frame{i}.f.^(1/2.2),['~/tbd_' int2str(i) frmt]);
	end
end
