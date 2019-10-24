mean_tiou = zeros(1,numel(seq));
mean_tiou_nc = zeros(1,numel(seq));
mean_tiou3d = zeros(1,numel(seq));
mean_tiou3d_nc = zeros(1,numel(seq));
mean_tiou3d_nc_oracle = zeros(1,numel(seq));
mean_tiou3d_nc3d_oracle = zeros(1,numel(seq));
mean_tiou_nc3d3d = zeros(1,numel(seq));
rerr = zeros(1,numel(seq));
rerr_gt = zeros(1,numel(seq));
rerr_est = zeros(1,numel(seq));
for i = 1:numel(seq)
	mean_tiou(i) = mean([tiou{i}(3:end)]);
	mean_tiou_nc(i) = mean([tiou_nc{i}(3:end)]);
	mean_tiou3d(i) = mean([tiou3d{i}(3:end)]);
	mean_tiou3d_nc(i) = mean([tiou3d_nc{i}(3:end)]);
	mean_tiou3d_nc_oracle(i) = mean([tiou3d_nc_oracle{i}(3:end)]);
	mean_tiou3d_nc3d_oracle(i) = nanmean([tiou3d_nc3d_oracle{i}(3:end)]);
	mean_tiou_nc3d3d(i) = mean([tiou_nc3d3d{i}(3:end)]);
	frms = [frame{i}{:}]; r_est = max([frms.Radius]);
	rr = [PAR{i}.R]; 
	if ~iscell(szs_gt{i})
		ind_r = linspace(1,numel(frame{i})+1,numel(rr));
		ind_use = linspace(3,numel(frame{i})+1,numel(rr));
		rerr(i) = mean(abs([rr(:)] - r_est));
		rerr_gt(i) = mean(abs(interp1(ind_r, [rr(:)], ind_use,'linear', 'extrap') - interp1(ind_gt{i}, szs_gt{i}, ind_use,'linear', 'extrap')));
		rerr_est(i) = mean(abs(interp1(ind_r, [rr(:)], ind_use,'linear', 'extrap') - interp1(ind_est{i}, szs_est{i}, ind_use,'linear', 'extrap')));
	else
		disp('Szs is cell');
	end
end

fprintf('[TbD] Mean TIoU  %.4f\n', mean(mean_tiou));
fprintf('[TbD-NC] Mean TIoU %.4f\n', mean(mean_tiou_nc));
fprintf('[TbD] Mean TIoU-3D %.4f\n', mean(mean_tiou3d));
fprintf('[TbD-NC] Mean TIoU-3D  %.4f\n', mean(mean_tiou3d_nc));
fprintf('[TbD-NC-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc_oracle));
fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc3d_oracle));
fprintf('[TbD-NC-3D] Mean TIoU-3D %.4f\n', mean(mean_tiou_nc3d3d));

fprintf('[TbD-NC] Mean absolute difference to GT radius %.4f\n', mean(rerr));
fprintf('[TbD-NC-3D] Mean absolute difference to GT radius %.4f\n', mean(rerr_est));
fprintf('[TbD-NC-3D-Oracle] Mean absolute difference to GT radius %.4f\n', mean(rerr_gt));


if false 
	rr = [PAR{i}.R]; plot(linspace(1,numel(frame{i}),numel(rr)), [rr(:)], 'g'); hold on; plot(ind_gt{i}, szs_gt{i}, 'r'); plot(ind_est{i}, szs_est{i}, 'b');
	plot(ind_gt{i},repmat(r_est,numel(ind_gt{i})),'m');
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));

	fprintf('[TbD-NC] Mean absolute difference to GT radius %s is %.4f\n', seq(i).name, rerr(i));
	fprintf('[TbD-NC-3D] Mean absolute difference to GT radiusfor %s is %.4f\n', seq(i).name, rerr_3d(i));
	fprintf('[TbD-NC-3D-Oracle] Mean absolute difference to GT radius for %s is %.4f\n', rerr_gt(i));
end