
params_tbd3d = [];
params_tbd3d.do_hier = true;
params_tbd3d.f0_maxiter = 10;
params_tbd3d.maxiter = 5;
params_tbd3d.iter_smoothing = 4;
params_tbd3d.do_intervals = true;
i = 1;

make_joint = true;

if false
	params_tbd3d.lambda_R = 0;
	params_tbd3d.alpha_cross_f = 0;
	params_tbd3d.alpha_cross_m = 0;
	params_tbd3d.do_hier = false;
end

[seq, folder] = EVAL.get_seq(0,'TbD-3D-aero');
n = 8;
resz = 1;

load('~/projects/data/TbD-3D-aero-norot-nocross-nohier.mat');

matF = matF_WB(matF_gt{1});
matM = matM_gt{1};

t = load(fullfile(folder, seq(i).name));
[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

bgr = median(Vk{i}, 4);

if make_joint
	iv = 1:size(Vk{i},4);
	ind = [14 20 37];
else
	iv = [12:14];
	ind_picked = [2 10 15 21];
	ind_picked_hs = [1 11 16 22];
end

% max_ind = max(iv)+1;
% [szs_gt{i},matF_gt{i},matM_gt{i},ind_gt{i}] = TD.estimate_3dtraj(im2double(Vk{i}(:,:,:,1:max_ind)), gt_coeffs{i}(1:max_ind), Ms{i}, n, params_tbd3d);

PAR_mod = PAR{i};
for jj = 1:numel(PAR_mod), PAR_mod(jj).R = PAR_mod(jj).R*1.1; end
[matF_hs, matM_hs] = TD.get_views_hs_3d(im2double(V_WB{i}),[],PAR_mod,n,true,false);

frmc = (iv(1) + iv(end))/2;


matF_inl = matF(:,:,:,ind_gt{i} >= iv(1) & ind_gt{i} <= iv(end)+1);
matM_inl = matM(:,:,:,ind_gt{i} >= iv(1) & ind_gt{i} <= iv(end)+1);
% matF_inl = matF_inl.*matM_inl;
matM_inl = matM_inl.^2.2;
matF_inl = matF_inl.^1.3;

rr = [PAR{i}.R]; 
ind_r = linspace(1,size(Vk{i},4)+1,numel(rr));

if make_joint
	for k = ind
		inp = [];
		f1 = matF_inl(:,:,:,ind_gt{i} >= k & ind_gt{i} < k+1);
		while size(f1,4) < 8
			f1 = f1(:,:,:,repelem(1:size(f1,4),2));
		end
		f1 = montage(f1,'Size',[1 size(f1,4)]); f1 = f1.CData;
		m1 = matM_inl(:,:,:,ind_gt{i} >= k & ind_gt{i} < k+1);
		imwrite(m1(:,:,:,1),['~/tmp/aero_Mone' int2str(k) '.png']);
		m1 = montage(m1,'Size',[1 size(m1,4)]); m1 = m1.CData;
		imwrite(f1,['~/tmp/aero_F' int2str(k) '.png']);
		imwrite(m1,['~/tmp/aero_M' int2str(k) '.png']);

		fgt = matF_hs(:,:,:,ind_r >= k & ind_r < k+1);
		fgt = montage(fgt,'Size',[1 size(fgt,4)]); fgt = fgt.CData;
		imwrite(fgt,['~/tmp/aero_FGT' int2str(k) '.png']);

		img = Vk{1}(:,:,:,k);
		inp = Vk_WB{1}(:,:,:,k);
		% new_coeff = [gt_coeffs{1}{k}{1}(:,1) sum(gt_coeffs{1}{k}{end}')'];
		new_coeff = [gt_coeffs{1}{k}{1}(:,1) gt_coeffs{1}{k}{end}(:,1)];
		new_coeff = {[new_coeff(:,1) new_coeff(:,2)-new_coeff(:,1)]};
		Hall = myTrajRender(size2(inp), new_coeff, [0 1]);
		inp(logical(cat(3,zeros(size(Hall)),Hall>0))) = 255;
		HM = conv2(Hall,ones(size(matM_inl(:,:,:,1))),'same');
		[yy,xx]=find(HM>0);
		inp = inp(min(yy(:)):max(yy(:)),min(xx(:)):max(xx(:)),:);
		imwrite(inp,['~/tmp/aero_in' int2str(k) '.png']);

		Hs = [];
		gap = 0.07;
		for k3 = 1:n
			H = evaluate_vcoeff(size2(img), gt_coeffs{1}{k}, [( ( (k3-1+gap) )/n) (k3/n)]);
			H = double(H);
			H = H ./ sum(H(:));
			Hs = cat(3, Hs, H);
		end
		Hs = Hs / sum(Hs(:));
		save(['~/tmp/aero_' int2str(k) '.mat'],'img','bgr','Hs','Hall','f1','m1','fgt');
	end	
else
	
	matF_hs_inl = matF_hs(:,:,:,ind_r >= iv(1) & ind_r <= iv(end)+1);

	matF_inl = matF_inl(:,:,:,ind_picked);
	matM_inl = matM_inl(:,:,:,ind_picked);
	matF_hs_inl = matF_hs_inl(:,:,:,ind_picked_hs);
	for frmi = iv
		imwrite(Vk_WB{i}(:,:,:,frmi),['~/tmp/aero_' int2str(frmi) '.png']);
	end

	for k = 1:numel(ind_picked)
		imwrite(matF_inl(:,:,:,k),['~/tmp/matF_' int2str(k) '.png']);
		imwrite(matM_inl(:,:,:,k),['~/tmp/matM_' int2str(k) '.png']);
		imwrite(matF_hs_inl(:,:,:,k),['~/tmp/matF_hs_' int2str(k) '.png']);
	end	
end