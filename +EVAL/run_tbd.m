if ~exist('params', 'var')
	[params, cfg] = EVAL.get_params(true);
end
[seq, folder] = EVAL.get_seq(0,'TbD');

if cfg.extend_dataset
	for i = 1:numel(seq)
		seq(i).end_ind = []; 
		seq(i).start_ind = 1; 
	end
end

tiou = repmat({[]}, 1, numel(seq));
tiou_nc = repmat({[]}, 1, numel(seq));
frame = repmat({[]}, 1, numel(seq));
curves = repmat({[]}, 1, numel(seq));
par = repmat({[]}, 1, numel(seq));

r0_all = double(zeros(1,numel(seq)));
exectime = zeros(size(seq));

% if isempty(gcp('nocreate')), parpool(numel(seq)); end
% parfor i = 1:numel(seq)
for i = 1:numel(seq)
 	warning('off','all');
	disp(seq(i).name);

	params0 = params;

	[~,sname,~] = fileparts(seq(i).name);
	sname = erase(sname,'-gc');
	tdir = load(fullfile(folder, 'templates', [sname '_template.mat']));
	if tdir.scale ~= seq(i).resize
		tdir.template = imresize(tdir.template, seq(i).resize / tdir.scale);
	end
	if params0.use_template	
		if ~isempty(seq(i).pad_pcent)
			[template m0] = setupTemplate(tdir.template,seq(i).pad_pcent);
		else
			[template m0] = setupTemplate(tdir.template);
		end
		params0.F = im2double(template).^(1/params0.gm);
		params0.M = m0;
	end
	r0_all(i) = double(size(tdir.template,1)/2);

	t = load(fullfile(folder, seq(i).name));

	if ~isempty(seq(i).end_ind)
		video = VID.VideoMat(imresize(t.Vk(:,:,:,seq(i).start_ind:seq(i).end_ind), seq(i).resize));
		par{i} = cellfun(@(x) x*seq(i).resize, t.PAR(seq(i).start_ind:seq(i).end_ind), 'UniformOutput', false);
	else
		video = VID.VideoMat(imresize(t.Vk(:,:,:,seq(i).start_ind:end), seq(i).resize));
		par{i} = cellfun(@(x) x*seq(i).resize, t.PAR(seq(i).start_ind:end), 'UniformOutput', false);
	end

	if params.apply_normalisation
		gm_val = 1;
		[counts,binLocations] = imhist(rgb2gray(mean(im2double(video.mat).^gm_val,4)));
		perc = counts / sum(counts);
		aperc = cumsum(perc);
		[~,ind] = find(aperc' < 0.999); ind = ind(end);
		val = max(0.5, binLocations(ind));
		video.mat = video.mat / val;
		params0.F = (params0.F.^gm_val) / val;
	end

	timev = tic;
	frame{i} = EVAL.tbd_main_loop(video, cfg, params0);
	exectime(i) = toc(timev);
	fprintf('Sequence %s took %.3f sec.\n',seq(i).short, exectime(i));  

	frame{i} = [repmat({Frame.empty},1,0) frame{i}];

	VIS.show_all(t, seq(i), frame{i});

	[tiou{i}] = FIT.gt_cost_iou(frame{i}, par{i}, r0_all(i));
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	[expf, ~] = estimate_exposure(frame{i});
	%% TbD-NC
	[curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6);
	[tiou_nc{i}] = FIT.gt_cost_iou_curves(curves{i}, frame{i}, par{i}, r0_all(i));
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
end


%% recalculate TIoU
for i = 1:numel(seq)
	if ~isempty(frame{i})
		[tiou{i}] = FIT.gt_cost_iou(frame{i}, par{i}, r0_all(i)); 
		[tiou_nc{i}] = FIT.gt_cost_iou_curves(curves{i}, frame{i}, par{i}, r0_all(i)); 
	end
end	

tiou_mean = cellfun(@(x) mean(x(3:end)), tiou);
tiou_nc_mean = cellfun(@(x) mean(x(3:end)), tiou_nc);
fprintf('TbD: '); fprintf('%.3f ', tiou_mean); fprintf('\n');
fprintf('TbD-NC: '); fprintf('%.3f ', tiou_nc_mean); fprintf('\n');
fprintf('[TbD] Mean TIoU is %.3f\n', mean(tiou_mean));
fprintf('[TbD-NC] Mean TIoU is %.3f\n', mean(tiou_nc_mean));

lens = cellfun(@(x) numel(x), frame);
meanex = exectime ./ lens; 	

