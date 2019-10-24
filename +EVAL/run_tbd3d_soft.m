[seq, folder] = EVAL.get_seq(0,'TbD-3D');
n = 8;
resz = 1;

if any(cellfun(@(x) isempty(x), matF))
	load('~/projects/data/TbD-3D.mat');
end

tiou = repmat({[]}, 1, numel(seq));
tiou_nc = repmat({[]}, 1, numel(seq));
tiou3d = repmat({[]}, 1, numel(seq));
tiou3d_nc = repmat({[]}, 1, numel(seq));
tiou3d_nc_oracle = repmat({[]}, 1, numel(seq));
tiou3d_nc3d_oracle = repmat({[]}, 1, numel(seq));
tiou_nc3d3d = repmat({[]}, 1, numel(seq));

matF_hs = repmat({[]}, 1, numel(seq));
matM_hs = repmat({[]}, 1, numel(seq));

for i = 1:numel(seq)
	disp(seq(i).name);

	if isempty(Vk_WB{i}) 
		t = load(fullfile(folder, seq(i).name));
		[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
		Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
		Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);
	end

	frms = [frame{i}{:}];
	r_est = max([frms.Radius]);

	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));

	%% TIoU-3D
	[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i}, 2);
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));

	[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2);
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));

	%% TIoU-3D TbD-Oracle
	[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));

	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind_gt{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));

	%% TbD-3D
	if isempty(curves{i})
		tiou_nc3d3d{i} = zeros(1,numel(frame{i}));
	else
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs_est{i}, ind_est{i});
	end
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));

	[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);
end

EVAL.get_stats;

if false
	PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i});
	PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i}, szs_gt{i}, ind_gt{i});
	PAPER.gen_tex_image_gt(curves{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i}, szs_est{i}, ind_est{i});
	% PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, szs_est{i}, ind_est{i});
	PAPER.generate_imgs(matF_gt{i}, frame{i}, matF_hs{i}, ind_gt{i}(1:8));
end
