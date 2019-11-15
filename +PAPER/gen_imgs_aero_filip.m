function [] = gen_imgs_aero_filip()
ind = [14 20 37];
n = 8;
resz = 1;
i = 1;
[seq, folder] = EVAL.get_seq(0,'TbD-3D-aero');
load('~/projects/data/TbD-3D-aero-norot-nocross-nohier.mat');
t = load(fullfile(folder, seq(i).name));
[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

rr = [PAR{i}.R]; 
ind_r = linspace(1,size(Vk{i},4)+1,numel(rr));
	
for k = ind
	load(['~/projects/data/aero_' int2str(k) '_res.mat']);
	mr = permute(mr,[1 2 4 3]);
	fr = permute(fr,[1 2 4 3]);

	matF = fr;
	matF(matF < 0) = 0; 
	matF(matF > 1) = 1;

	matF = matF_WB(fr);
	matF = matF.^1.3;
	
	% matF = matF.^0.4;

	matM = mr;
	% matM = mr.^1.3;

	f1 = montage(matF,'Size',[1 size(matF,4)]); f1 = f1.CData;
	m1 = montage(matM,'Size',[1 size(matM,4)]); m1 = m1.CData;
	imwrite(f1,['~/tmp/aerof_F' int2str(k) '.png']);
	imwrite(m1,['~/tmp/aerof_M' int2str(k) '.png']);
	
	PAR_mod = PAR{i};
	for jj = 1:numel(PAR_mod), PAR_mod(jj).R = repmat(size(matF,1)/2, [1 numel(PAR_mod(jj).R)]); end
	[matF_hs, matM_hs] = TD.get_views_hs_3d(im2double(V_WB{i}),[],PAR_mod,n,true,false);

	fgt = matF_hs(:,:,:,ind_r >= k & ind_r < k+1);
	fgt = montage(fgt,'Size',[1 size(fgt,4)]); fgt = fgt.CData;
	imwrite(fgt,['~/tmp/aerof_FGT' int2str(k) '.png']);

end