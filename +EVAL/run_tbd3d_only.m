
test_oracle = false;
test_tbd3d = false;

if isempty(params_tbd3d)
	params_tbd3d.do_hier = true;
	params_tbd3d.f0_maxiter = 10;
	params_tbd3d.maxiter = 5;
	params_tbd3d.iter_smoothing = 4;
end

% if isempty(gcp('nocreate')), parpool(numel(seq)); end
% parfor i = 1:numel(seq)
for i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name);

	if test_oracle
		[szs_gt{i},matF_gt{i},matM_gt{i},ind_gt{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
		[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind_gt{i});
		fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));
	end

	if test_tbd3d
		[szs_est{i},matF{i},matM{i},ind_est{i}] = TD.estimate_3dtraj(im2double(Vk{i}), curves{i}, Ms{i}, n, params_tbd3d);
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs_est{i}, ind_est{i});
		fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
	end
end

EVAL.get_stats;

if false
	save('~/projects/data/TbD-3D.mat','frame','curves','Ms','gt_coeffs','ind_est','ind_gt','matF_gt','matM_gt','szs_gt','szs_est','matF','matM');
	PAPER.print_table_overleaf(PAPER.generate_table_comparison());
end
