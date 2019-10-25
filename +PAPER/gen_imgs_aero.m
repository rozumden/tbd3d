
params_tbd3d = [];
params_tbd3d.do_hier = true;
params_tbd3d.f0_maxiter = 10;
params_tbd3d.maxiter = 5;
params_tbd3d.iter_smoothing = 4;
params_tbd3d.do_intervals = true;
i = 1;

[seq, folder] = EVAL.get_seq(0,'TbD-3D-aero');
n = 8;
resz = 1;

load('~/projects/data/TbD-3D-aero.mat');

% [szs_gt{i},matF_gt{i},matM_gt{i},ind_gt{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);

matF = matF_WB(matF_gt{1});
matM = matM_gt{1};

t = load(fullfile(folder, seq(i).name));
[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

PAR_mod = PAR{i};
% for jj = 1:numel(PAR_mod), PAR_mod(jj).R = PAR_mod(jj).R*1.1; end
[matF_hs, matM_hs] = TD.get_views_hs_3d(im2double(V_WB{i}),[],PAR_mod,n,true);

iv = [12:14];
frmc = (iv(1) + iv(end))/2;

for frmi = iv
	imwrite(Vk_WB{i}(:,:,:,frmi),['aero_' int2str(frmi) '.png']);
end

matF_inl = matF(:,:,:,ind_gt{i} >= iv(1) & ind_gt{i} <= iv(end)+1);
matM_inl = matM(:,:,:,ind_gt{i} >= iv(1) & ind_gt{i} <= iv(end)+1);
matF_inl = matF_inl.*matM_inl;
matM_inl = matM_inl.^2.2;

ind_picked = [2 10 15 21];

matF_inl = matF_inl(:,:,:,ind_picked);
matM_inl = matM_inl(:,:,:,ind_picked);


rr = [PAR{i}.R]; 
ind_r = linspace(1,size(Vk{i},4)+1,numel(rr));
matF_hs_inl = matF_hs(:,:,:,ind_r >= iv(1) & ind_r <= iv(end)+1);

ind_picked_hs = [1 11 16 22];

matF_hs_inl = matF_hs_inl(:,:,:,ind_picked_hs);

for k = 1:numel(ind_picked)
	imwrite(matF_inl(:,:,:,k),['matF_' int2str(k) '.png']);
	imwrite(matM_inl(:,:,:,k),['matM_' int2str(k) '.png']);
	imwrite(matF_hs_inl(:,:,:,k),['matF_hs_' int2str(k) '.png']);
end

