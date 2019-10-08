
test_oracle = true;
test_tbd3d = false;

if isempty(params_tbd3d)
	params_tbd3d.do_hier = true;
	params_tbd3d.f0_maxiter = 10;
	params_tbd3d.maxiter = 5;
	params_tbd3d.iter_smoothing = 4;
end

if isempty(gcp('nocreate')), parpool(numel(seq)); end
parfor i = 1:numel(seq)
	warning('off','all');
	% disp(seq(i).name);

	if test_oracle
		[szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
		[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i});
		fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));
	end

	if test_tbd3d
		if isempty(curves{i})
			tiou_nc3d3d{i} = zeros(1,numel(frame{i}));
		else
			[szs{i},matF{i},matM{i}] = TD.estimate_3dtraj(im2double(Vk{i}), curves{i}, Ms{i}, n);
			[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs{i});
		end
		fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
	end
end


for i = 1:numel(seq)
	mean_tiou3d_nc_oracle(i) = mean([tiou3d_nc_oracle{i}(3:end)]);
	mean_tiou3d_nc3d_oracle(i) = mean([tiou3d_nc3d_oracle{i}(3:end)]);
	mean_tiou_nc3d3d(i) = mean([tiou_nc3d3d{i}(3:end)]);
	frms = [frame{i}{:}]; r_est = max([frms.Radius]);
	rr = [PAR{i}.R]; 
	ss_gt = [szs_gt{i}{:}]; 
	ss_3d = r_est;
	if ~isempty(szs{i}), ss_3d = [szs{i}{:}]; end
	rerr(i) = mean(abs([rr(:)] - r_est));
	rerr_gt(i) = mean(abs([rr(:)] - [ss_gt(:)]));
	rerr_3d(i) = mean(abs([rr(:)] - [ss_3d(:)]));
end

fprintf('[TbD-NC] Mean absolute difference to GT radius %.4f\n', mean(rerr));
fprintf('[TbD-NC-3D] Mean absolute difference to GT radius %.4f\n', mean(rerr_3d));
fprintf('[TbD-NC-3D-Oracle] Mean absolute difference to GT radius %.4f\n', mean(rerr_gt));


if test_oracle
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc_oracle));
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc3d_oracle));
end

if test_tbd3d
	fprintf('[TbD-NC] Mean TIoU-3D  %.4f\n', mean(mean_tiou3d_nc));
	fprintf('[TbD-NC-3D] Mean TIoU-3D %.4f\n', mean(mean_tiou_nc3d3d));
end

if false
	PAPER.print_table_overleaf(PAPER.generate_table_comparison());
end
