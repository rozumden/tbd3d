function [] = gen_cvpr20_supp_falling()
n = 10;
resz = 6;
i = 2;
k_before = 0;
k_after = 0;


[seq, folder] = EVAL.get_seq(0,'falling');

load(['~/projects/data/falling_res2.mat']);
t = load(fullfile(folder, seq(i).name([1:end-11 end-8:end])));
Vk{i} = permute(t.V(:,end:-1:1,:,:),[2 1 3 4]);
tgt = load(fullfile(folder, seq(i).name));
HS{i} = permute(tgt.HS(:,end:-1:1,:,:),[2 1 3 4]);

nimgs = 3;

for k = ki
	if nimgs == 2
		writer = VID.VideoWriterWrapperDoubleImg('~/projects/vis/', ['key' int2str(k)], 'fps', 10);
	elseif nimgs == 3
		writer = VID.VideoWriterWrapperTriple('~/projects/vis/', ['key' int2str(k)], 'fps', 10);
	else nimgs == 4
		writer = VID.VideoWriterWrapperQuadruple('~/projects/vis/', ['key' int2str(k)], 'fps', 20);
	end

	bgr = median(im2double(Vk{i}),4);
	Hall = evaluate_vcoeff(size2(bgr), gt_coeffs{i}{1}, [0 1]);

	HHH = conv2(Hall, ones(80),'same');
	[yy,xx] = find(HHH > 0);
	ys = [min(xx(:)) max(xx(:))];
	xs = [min(yy(:)) max(yy(:))];

	imin = Vk{i}(:,:,:,k);
	imin = imin(xs(1):xs(2),ys(1):ys(2),:);

	% kk1 = [1 2 4 5 6 7 9 10];
	kk1 = [1 2 2 3 4 5 6 7 7 8];
	for kk = 1:n
		in = n*(k-1) + kk;
		im_gt = HS{i}(:,:,:,in);

		H = evaluate_vcoeff(size2(bgr), gt_coeffs{i}{1}, [(kk-1)/n (kk/n)]);
		H = H / sum(H(:));
		
		F = matF_gt{i}(:,:,:,kk1(kk));
		M = matM_gt{i}(:,:,:,kk1(kk));

		F(F < 0) = 0; F(F > 1) = 1;
		M(M < 0) = 0; M(M > 1) = 1;

		F(~repmat(M > 0,[1 1 3])) = 0; 
		HF = cat(3, conv2(H, F(:,:,1), 'same'), conv2(H, F(:,:,2), 'same'),conv2(H, F(:,:,3), 'same'));
		HM = conv2(H, M, 'same');
		im = HF + (1 - HM).*bgr;

		HMvis = repmat(HM(xs(1):xs(2),ys(1):ys(2),:),[1 1 3]);
		HMvis = HMvis/max(HMvis(:));
		im_gt = im_gt(xs(1):xs(2),ys(1):ys(2),:);
		im = im(xs(1):xs(2),ys(1):ys(2),:);

		if nimgs == 4
			writer.write_buffer_img(imin);
			writer.write_buffer2_img(255*im);
			writer.write_buffer3_img(255*HMvis);
			writer.write_img(im_gt);
		elseif nimgs == 3
			writer.write_buffer_img(imresize(imin,resz))
			writer.write_buffer2_img(imresize(255*im,resz));
			writer.write_img(imresize(im_gt,resz));
		elseif nimgs == 2
			writer.write_buffer_img(imresize(imin,resz));
			writer.write_img(imresize(255*im,resz));
		end
	end
	writer.close();

end