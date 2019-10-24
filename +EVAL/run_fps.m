resz = 0.5;

params_tbd3d = [];
params_tbd3d.do_hier = true;
params_tbd3d.f0_maxiter = 10;
params_tbd3d.maxiter = 5;
params_tbd3d.iter_smoothing = 4;
params_tbd3d.do_intervals = true;

[params, cfg] = EVAL.get_params(true);
[seq, folder] = EVAL.get_seq(0,'TbD-3D');

averages = struct('tiou',[],'tiou_nc',[],'tiou3d',[],'tiou3d_nc',[],'tiou3d_nc_oracle',[],'tiou3d_nc3d_oracle',[],'nc3d3d',[],'rerr',[],'rerr_est',[],'rerr_gt',[]);

ns = [8 4 2 1];
for n = ns
	data = load(['~/projects/data/TbD-3D-n' int2str(n) '.mat']);
	frame = data.frame; curves = data.curves;
	parfor i = 1:numel(seq)
		disp(seq(i).name);
		% try
			t = load(fullfile(folder, seq(i).name));
			[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
			Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
			Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);

			frms = [data.frame{i}{:}];
			r_est = max([frms.Radius]);
			Ms{i} = double(diskMask([],r_est+5));

			[tiou{i}] = FIT.gt_cost_iou_3d(data.frame{i}, PAR{i});
			[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(data.curves{i}, data.frame{i}, PAR{i});

			[tiou3d{i}] = FIT.gt_cost_iou_3d(data.frame{i}, PAR{i}, 2);
			[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(data.curves{i}, data.frame{i}, PAR{i}, 2);
			[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});

			[szs_gt{i},matF_gt{i},matM_gt{i},ind_gt{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
			[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind_gt{i});

			[szs_est{i},matF{i},matM{i},ind_est{i}] = TD.estimate_3dtraj(im2double(Vk{i}), data.curves{i}, Ms{i}, n, params_tbd3d);
			[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(data.curves{i}, data.frame{i}, PAR{i}, 2, szs_est{i}, ind_est{i});
		% catch err
		% 	disp(['Error at '  seq(i).name]);
		% 	err
		% end
		% [matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);
		disp(['Finished ' seq(i).name]);
	end

	EVAL.get_stats;

	averages.tiou = [averages.tiou mean(mean_tiou)];
	averages.tiou_nc = [averages.tiou_nc mean(mean_tiou_nc)];
	averages.tiou3d = [averages.tiou3d mean(mean_tiou3d)];
	averages.tiou3d_nc = [averages.tiou3d_nc mean(mean_tiou3d_nc)];
	averages.tiou3d_nc_oracle = [averages.tiou3d_nc_oracle mean(mean_tiou3d_nc_oracle)];
	averages.tiou3d_nc3d_oracle = [averages.tiou3d_nc3d_oracle mean(mean_tiou3d_nc3d_oracle)];
	averages.nc3d3d = [averages.nc3d3d mean(mean_tiou_nc3d3d)];

	averages.rerr = [averages.rerr mean(rerr)];
	averages.rerr_est = [averages.rerr_est mean(rerr_est)];
	averages.rerr_gt = [averages.rerr_gt mean(rerr_gt)];

	save(['~/projects/data/TbD-3D-n' int2str(n) '_post.mat'],'averages','szs_gt','ind_gt','matF_gt','matM_gt','szs_est','ind_est','matF','matM');

end
