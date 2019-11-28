function [] = gen_cvpr20_supp_deformation()
ind = [14 20 37];
% ind = [14];
n = 8;
resz = 1;
i = 1;
k_before = 2;
k_after = 2;

[seq, folder] = EVAL.get_seq(0,'TbD-3D-aero');
load('~/projects/data/TbD-3D-aero-norot-nocross-nohier.mat');
t = load(fullfile(folder, seq(i).name));
[Vk{i},Vk_WB{i},PAR{i},V_WB{i},V_WB_POS] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

rr = [PAR{i}.R]; 
ind_r = linspace(1,size(Vk{i},4)+1,numel(rr));

for k = ind
	load(['~/projects/data/aero_' int2str(k) '_res.mat']);
	mr = permute(mr,[1 2 4 3]);
	fr = permute(fr,[1 2 4 3]);

	matF = fr;

	% matF(matF < 0) = 0; 
	% matF(matF > 1) = 1;
	% matF = matF_WB(fr);
	% matF = matF.^1.3;

	matM = mr;

	% keyboard
		
	writer = VID.VideoWriterWrapperQuadruple('~/projects/vis/', ['deformation' int2str(k) '_4videos'], 'fps', 20);


	bgr = median(im2double(Vk_WB{i}),4);
	Hall = zeros(size2(bgr));
	for ki = k-k_before:k+k_after
		Hall = Hall + evaluate_vcoeff(size2(bgr), gt_coeffs{i}{ki}, [0 1]);
	end
	HHH = conv2(Hall, ones(80),'same');
	[yy,xx] = find(HHH > 0);
	ys = [min(xx(:)) max(xx(:))];
	xs = [min(yy(:)) max(yy(:))];
	keyboard
	for ki = k-k_before:k+k_after
		imin = Vk_WB{i}(:,:,:,ki);
		imin = imin(xs(1):xs(2),ys(1):ys(2),:);
	
		for kk = 1:n
			in = n*(ki-1) + kk;
			im_gt = V_WB{i}(:,:,:,in);

			[~,in_gt] = min(abs(ind_gt{i} - ind_r(in)));

			H = evaluate_vcoeff(size2(bgr), gt_coeffs{i}{ki}(kk), [0 1]);
			H = H / sum(H(:));

			
			F = matF_gt{i}(:,:,:,in_gt);
			M = matM_gt{i}(:,:,:,in_gt);

			if k == ki
				% F = matF(:,:,:,kk);
				% M = matM(:,:,:,kk); 
			end
			% F = F.^0.5;

			F(F < 0) = 0; F(F > 1) = 1;
			M(M < 0) = 0; M(M > 1) = 1;
			F = matF_WB(F);
			F = F.^1.3;

			F(~repmat(M > 0,[1 1 3])) = 0; 
			HF = cat(3, conv2(H, F(:,:,1), 'same'), conv2(H, F(:,:,2), 'same'),conv2(H, F(:,:,3), 'same'));
			HM = conv2(H, M, 'same');
			im = HF + (1 - HM).*bgr;

			HMvis = repmat(HM(xs(1):xs(2),ys(1):ys(2),:),[1 1 3]);
			HMvis = HMvis/max(HMvis(:));
			im_gt = im_gt(xs(1):xs(2),ys(1):ys(2),:);
			im = im(xs(1):xs(2),ys(1):ys(2),:);

			writer.write_buffer_img(imin);
			writer.write_buffer2_img(255*im);
			writer.write_buffer3_img(255*HMvis);
			writer.write_img(im_gt);
		end
	end
	writer.close();

end