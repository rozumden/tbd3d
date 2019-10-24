function [] = generate_imgs(matF, frame, matF_hs, iv)
if nargin < 4
	iv = 41:0.5:43;
end
fc = [1.9 1 1.8];
WB = [2 1 2]; gamma_coef = 0.4;
for k = 1:3
	matF(:,:,k,:) = matF(:,:,k,:) ./ fc(k);
end

n = size(matF_hs,4) / numel(frame);
frmt = '.png';

for i = iv
	it = (i-1)*n + 1;
	F3D = matF(:,:,:,it); 
	F3D = ((F3D.*reshape(WB,1,1,[])/(max(WB))).^gamma_coef);

	GT = matF_hs(:,:,:,(i-1)*n + 1); 
	imwrite(GT,['~/gt_' int2str(it) frmt]);
	imwrite(F3D,['~/tbd3d_' int2str(it) frmt]);
	if round(i) == i
		imwrite(frame{i}.f.^(1/2.2),['~/tbd_' int2str(it) frmt]);
	end
end
