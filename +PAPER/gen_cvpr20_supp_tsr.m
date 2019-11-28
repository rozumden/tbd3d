n = 8;
resz = 1;
i = 1;
% i = 7;
writer = VID.VideoWriterWrapperTriple('~/projects/vis/', ['tsr' int2str(i)], 'fps', 20);

% PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}));

rr = [PAR{i}.R]; 
ind_r = linspace(1,size(Vk{i},4)+1,numel(rr)+1); ind_r = ind_r(1:(end-1));
bgr = median(im2double(Vk_WB{i}),4);
for in = 1:numel(ind_r)
	im_gt = V_WB{i}(:,:,:,in);
	k = floor(ind_r(in));
	kk = n*(ind_r(in) - k) + 1;

	[~,in_gt] = min(abs(ind_gt{i} - ind_r(in)));

	H = evaluate_vcoeff(size2(bgr), gt_coeffs{i}{k}(kk), [0 1]);
	H = H / sum(H(:));
	F = matF_gt{i}(:,:,:,in_gt);
	F(F < 0) = 0; F(F > 1) = 1;
	F = F.^0.5;
	M = matM_hs{i}(:,:,:,in); 
	M = imresize(M, size2(matM_gt{i}));
	M(M < 0) = 0; M(M > 1) = 1;
	F(~repmat(M > 0,[1 1 3])) = 0; 
	HF = cat(3, conv2(H, F(:,:,1), 'same'), conv2(H, F(:,:,2), 'same'),conv2(H, F(:,:,3), 'same'));
	HM = conv2(H, M, 'same');
	im = HF + (1 - HM).*bgr;

	imin = Vk_WB{i}(:,:,:,k);
	writer.write_buffer_img(imin);
	writer.write_buffer2_img(255*im);
	writer.write_img(im_gt);

end

writer.close();