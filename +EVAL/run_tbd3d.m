if ~exist('params', 'var')
	[params, cfg] = EVAL.get_params(true);
end
params.use_template = true;
params.alpha_F = 0.5;
[seq, folder] = EVAL.get_seq(0,'TbD-3D');
k = 8;
resz = 0.5;

for i = 4
	disp(seq(i).name);

	t = load(fullfile(folder, seq(i).name));
	[~,Vk,PAR] = generate_lowFPSvideo(t.V,t.POS,t.R,k,resz);

	params.BGR = median(im2double(Vk),4);
	params.M = double(diskMask([],max(max([PAR.R]))));
	params.F = [];
	video = VID.VideoMat(Vk);

	timev = tic;
	frame{i} = EVAL.tbd_main_loop(video, cfg, params);
	exectime(i) = toc(timev);
	fprintf('Sequence %s took %.3f sec.\n',seq(i).name, exectime(i));  

	VIS.show_all([], seq(i), frame{i});

	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR);
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	% [expf, ~] = estimate_exposure(frame{i}); 
	expf = 1;
	
	%% TbD-NC
	[curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6);
	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR);
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));

	[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR, 2);
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));

	[tiou_nc3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR, 2);
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d{i}(3:end)]));

	%% TbD-3D
	frms = [frame{i}{:}];
	r = max([frms.Radius]);
	M = double(diskMask([],r));
	[szs,matF,matM] = TD.estimate_3dtraj(im2double(Vk), curves{i}, M);

	[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR, 2);
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));
end


