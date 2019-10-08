if ~exist('params', 'var')
	[params, cfg] = EVAL.get_params(true);
end
params.use_template = true;
params.alpha_F = 1;
[seq, folder] = EVAL.get_seq(0,'TbD-3D');
n = 8;
resz = 0.5;

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

if ~exist('Vk','var')
	Vk = repmat({[]}, 1, numel(seq));
	PAR = repmat({[]}, 1, numel(seq));
	V_WB = repmat({[]}, 1, numel(seq));
	Ms = repmat({[]}, 1, numel(seq));
end

matF_hs = repmat({[]}, 1, numel(seq));
matM_hs = repmat({[]}, 1, numel(seq));
matF_gt = repmat({[]}, 1, numel(seq));
matM_gt = repmat({[]}, 1, numel(seq));
matF = repmat({[]}, 1, numel(seq));
matM = repmat({[]}, 1, numel(seq));
ind = repmat({[]}, 1, numel(seq));
szs = repmat({[]}, 1, numel(seq));
szs_gt = repmat({[]}, 1, numel(seq));

for i = 4
% for i = 1:numel(seq)
% if isempty(gcp('nocreate')), parpool(numel(seq)); end
% parfor i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name);
	params0 = params;

	if isempty(Vk{i})
		t = load(fullfile(folder, seq(i).name));
		[~,Vk{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
	end

	params0.BGR = median(im2double(Vk{i}),4);
	params0.M = double(diskMask([],max(max([PAR{i}.R]))));
	params0.F = [];
	video = VID.VideoMat(Vk{i});

	timev = tic;
	frame{i} = EVAL.tbd_main_loop(video, cfg, params0);
	exectime(i) = toc(timev);
	fprintf('Sequence %s took %.3f sec.\n',seq(i).name, exectime(i));  

	VIS.show_all([], seq(i), frame{i});

	% [expf, ~] = estimate_exposure(frame{i}); 
	expf = 1;
	frms = [frame{i}{:}];
	r_est = max([frms.Radius]);
	Ms{i} = double(diskMask([],r_est));

	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	[curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6);
	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));

	% VIS.curve(curves{i});

	%% TIoU-3D
	[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i}, 2);
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));

	[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2);
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));

	%% TIoU-3D TbD-Oracle
	[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));

	[szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n);
	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));

	%% TbD-3D
	if isempty(curves{i})
		tiou_nc3d3d{i} = zeros(1,numel(frame{i}));
	else
		[szs{i},matF{i},matM{i}] = TD.estimate_3dtraj(im2double(Vk{i}), curves{i}, Ms{i}, n);
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs{i});
	end
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));

	[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},Ms{i},n,true);
end
	% catch
	% 	disp(['Error in ' int2str(i) ' : ' seq(i).name]);
	% end
if false
	save('~/projects/data/TbD-3D.mat','Vk','PAR','V_WB','frame','curves','Ms','gt_coeffs','ind','matF_gt','matM_gt','szs_gt','szs','matF','matM');
end


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


if false 
	rr = [PAR{i}.R]; plot(ind{i}, [rr(:)], 'g'); hold on; ss1 = [szs_gt{i}{:}]; plot(ind{i}, [ss1(:)], 'r'); ss2 = [szs{i}{:}]; plot(ind{i}, [ss2(:)], 'b');
	plot(ind{i},repmat(r_est,numel(ind{i})),'m');
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
end