function [] = gen_imgs_tennis_filip()
ind = 58;
n = 8;
resz = 1;
i = 4;
[seq, folder] = EVAL.get_seq(0,'TbD');
load('../data/TbD-T1(0.713)-LT.mat')

t = load(fullfile(folder, seq(i).name));
[gt_coeffs{i}] = par2coeff(t.PAR);
	
for id = ind
	load(['~/projects/data/tennis_res.mat']);
	mr = permute(mr,[1 2 4 3]);
	fr = permute(fr,[1 2 4 3]);

	matF = fr;
	% matF(matF < 0) = 0; 
	% matF(matF > 1) = 1;

	% matF = matF_WB(fr);
	% matF = matF.^1.3;
	
	% matF = matF.^0.4;

	matM = mr;
	% matM = mr.^1.3;

	f1 = montage(matF,'Size',[1 size(matF,4)]); f1 = f1.CData;
	m1 = montage(matM,'Size',[1 size(matM,4)]); m1 = m1.CData;
	imwrite(f1,['~/tmp/tennis_F.png']);
	imwrite(m1,['~/tmp/tennis_M.png']);
	
	r = size(f1,1)/2;
	fgt = [];
	for k = 1:8
		ind = (id-1)*n + k + 1;
		nma = sprintf('/mnt/lascar/rozumden/dataset/TbD_hs/tennis/%08d.png', ind);
		img = imread(nma);
		ps = gt_coeffs{i}{id}{k}(:,1) + gt_coeffs{i}{id}{k}(:,2);
		f = img(ps(2)-r:ps(2)+r,ps(1)-r:ps(1)+r,:);
		fgt = cat(2,fgt,f);
	end
	imwrite(fgt,['~/tmp/tennis_FGT.png']);
end