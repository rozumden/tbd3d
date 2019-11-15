[params, cfg] = EVAL.get_params(true);
[params_tbd3d] = EVAL.get_3dparams();
[seq, folder] = EVAL.get_seq(0,'gym');
n = 8;
resz = 0.5;

for i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name);

	t = load(fullfile(folder, seq(i).name));
	[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
	Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
	Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

	r_est = max(max([PAR{i}.R]));
	Ms{i} = double(diskMask([],r_est));

	[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));

	[szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));

	[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);

end

EVAL.run_tbd3d_rotation;


for i = 1:numel(seq)
	disp(seq(i).name);
	[~,sname,~] = fileparts(seq(i).name);
	sname_gt = [sname '_Phys.mat'];
	t_gt = load(fullfile(folder, sname_gt));
	
	t_gt.Omega_GT;
	t_gt.p{3};

end

EVAL.get_stats;
