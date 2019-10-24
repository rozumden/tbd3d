[params, cfg] = EVAL.get_params(true);

params_tbd3d = [];
params_tbd3d.do_hier = true;
params_tbd3d.f0_maxiter = 10;
params_tbd3d.maxiter = 5;
params_tbd3d.iter_smoothing = 4;
params_tbd3d.do_intervals = true;

[seq, folder] = EVAL.get_seq(0,'TbD-3D-aero');
n = 8;
resz = 1;

frame = repmat({[]}, 1, numel(seq));
curves = repmat({[]}, 1, numel(seq));
exectime = zeros(size(seq));
gt_coeffs = repmat({[]}, 1, numel(seq));
tiou = repmat({[]}, 1, numel(seq));
tiou_nc = repmat({[]}, 1, numel(seq));
tiou3d = repmat({[]}, 1, numel(seq));
tiou3d_nc = repmat({[]}, 1, numel(seq));
tiou3d_nc_oracle = repmat({[]}, 1, numel(seq));
tiou3d_nc3d_oracle = repmat({[]}, 1, numel(seq));
tiou_nc3d3d = repmat({[]}, 1, numel(seq));

Vk = repmat({[]}, 1, numel(seq));
Vk_WB = repmat({[]}, 1, numel(seq));
PAR = repmat({[]}, 1, numel(seq));
V_WB = repmat({[]}, 1, numel(seq));
Ms = repmat({[]}, 1, numel(seq));

matF_hs = repmat({[]}, 1, numel(seq));
matM_hs = repmat({[]}, 1, numel(seq));
matF_gt = repmat({[]}, 1, numel(seq));
matM_gt = repmat({[]}, 1, numel(seq));
matF = repmat({[]}, 1, numel(seq));
matM = repmat({[]}, 1, numel(seq));
ind_gt = repmat({[]}, 1, numel(seq));
ind_est = repmat({[]}, 1, numel(seq));
szs_gt = repmat({[]}, 1, numel(seq));
szs_est = repmat({[]}, 1, numel(seq));

do_tbd = false;

for i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name);
	params0 = params;

	if isempty(Vk{i})
		t = load(fullfile(folder, seq(i).name));
		[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
		Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
		Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);
	end

	params0.BGR = median(im2double(Vk{i}),4);
	params0.M = double(diskMask([],max(max([PAR{i}.R]))));
	params0.F = repmat(params0.M,[1 1 3]);
	video = VID.VideoMat(Vk{i});

	Ms{i} = double(diskMask([],max(max([PAR{i}.R]))));

	if do_tbd
		timev = tic;
		frame{i} = EVAL.tbd_main_loop(video, cfg, params0);
		exectime(i) = toc(timev);
		fprintf('Sequence %s took %.3f sec.\n',seq(i).name, exectime(i));  

		VIS.show_all([], seq(i), frame{i});

		% [expf, ~] = estimate_exposure(frame{i}); 
		expf = 1;
		frms = [frame{i}{:}];
		r_est = max([frms.Radius]);
		Ms{i} = double(diskMask([],r_est+5));

		[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
		fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

		[curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6);
		[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
		fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
	end

	% VIS.curve(curves{i});
	% VIS.curve(gt_coeffs{i});

	%% TIoU-3D
	if do_tbd
		[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i}, 2);
		fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));

		[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2);
		fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));
	end

	%% TIoU-3D TbD-Oracle
	[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));

	[szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));

	%% TbD-3D
	if do_tbd
		[szs_est{i},matF{i},matM{i},ind_est{i}] = TD.estimate_3dtraj(im2double(Vk{i}), curves{i}, Ms{i}, n, params_tbd3d);
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs_est{i}, ind_est{i});
		fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
	end

	[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);
end

if false
	save('~/projects/data/TbD-3D-aero.mat','frame','curves','Ms','gt_coeffs','ind_est','ind_gt','matF_gt','matM_gt','szs_gt','szs_est','matF','matM');
end

EVAL.get_stats;
