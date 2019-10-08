function [srctable] = generate_table_comparison()
load('~/projects/data/TbD-3D.mat');

precision = 3;

sep = ' & ';
srctable = '';

[seq, folder] = EVAL.get_seq(0,'TbD-3D');
for i = 1:numel(seq)
	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));
	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
	[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i}, 2);
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));
	[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2);
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));
	frms = [frame{i}{:}]; r_est = max([frms.Radius]);
	[tiou3d_nc_oracle{i},~,~] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));
	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));
	if isempty(curves{i})
		tiou_nc3d3d{i} = zeros(1,numel(frame{i}));
	else
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs{i});
	end
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
	% [matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},Ms{i},n,true);
end

disp('--------------------------------------------------');

mean_tiou = zeros(1,numel(seq));
mean_tiou_nc = zeros(1,numel(seq));
mean_tiou3d = zeros(1,numel(seq));
mean_tiou3d_nc = zeros(1,numel(seq));
mean_tiou3d_nc_oracle = zeros(1,numel(seq));
mean_tiou3d_nc3d_oracle = zeros(1,numel(seq));
mean_tiou_nc3d3d = zeros(1,numel(seq));
rerr = zeros(1,numel(seq));
rerr_gt = zeros(1,numel(seq));
rerr_3d = zeros(1,numel(seq));
for i = 1:numel(seq)
	mean_tiou(i) = mean([tiou{i}(3:end)]);
	mean_tiou_nc(i) = mean([tiou_nc{i}(3:end)]);
	mean_tiou3d(i) = mean([tiou3d{i}(3:end)]);
	mean_tiou3d_nc(i) = mean([tiou3d_nc{i}(3:end)]);
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

fprintf('[TbD] Mean TIoU  %.4f\n', mean(mean_tiou));
fprintf('[TbD-NC] Mean TIoU %.4f\n', mean(mean_tiou_nc));
fprintf('[TbD] Mean TIoU-3D %.4f\n', mean(mean_tiou3d));
fprintf('[TbD-NC] Mean TIoU-3D  %.4f\n', mean(mean_tiou3d_nc));
fprintf('[TbD-NC-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc_oracle));
fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D %.4f\n', mean(mean_tiou3d_nc3d_oracle));
fprintf('[TbD-NC-3D] Mean TIoU-3D %.4f\n', mean(mean_tiou_nc3d3d));

fprintf('[TbD-NC] Mean absolute difference to GT radius %.4f\n', mean(rerr));
fprintf('[TbD-NC-3D] Mean absolute difference to GT radius %.4f\n', mean(rerr_3d));
fprintf('[TbD-NC-3D-Oracle] Mean absolute difference to GT radius %.4f\n', mean(rerr_gt));

tious = [mean_tiou; mean_tiou_nc; mean_tiou3d; mean_tiou3d_nc; mean_tiou_nc3d3d; mean_tiou3d_nc_oracle; mean_tiou3d_nc3d_oracle];
names = {'TbD', 'TbD-NC', 'TbD', 'TbD-NC', 'TbD-3D', 'TbD-Oracle', 'TbD-3D-Oracle'};

for i = 1:numel(seq)
	[~,seqname,~] = fileparts(seq(i).name);
	C = strsplit(seqname,'_');
	short = C{end};

	c3 = int2str(numel(frame{i}));
	srctable = [srctable short sep c3];

	[~,ind] = max(tious(:,i)); ind = find(tious(ind,i) == tious(:,i));
	for k = 1:numel(names)
		col = sprintf(['%.' int2str(precision) 'f'],tious(k,i));

		lf = ''; rg = '';
		if any(k == ind)
			lf = '\\textbf{';
			rg = '}';
		end

		srctable = [srctable sep lf col(2:end) rg];
	end
	srctable = [srctable '\\\\  \n '];
end

meanfrm = sprintf('%.0f',mean(arrayfun(@(x) numel(x{1}), frame)));
srctable = [srctable '\\hline \n  Average ' sep meanfrm];

[~,ind] = max(mean(tious'));
for k = 1:numel(names)
	meaniou = sprintf(['%.' int2str(precision) 'f'],mean(tious(k,:)));
	lf = ''; rg = '';
	if k == ind
		lf = '\\textbf{';
		rg = '}';
	end

	srctable = [srctable sep lf meaniou(2:end) rg];
end

srctable = [srctable '\\\\\n'];

