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
for i = 4
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

	frame{i} = [repmat({VID.Frame.empty},1,0) frame{i}];

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



if false
	[gt_coeffs{i}] = par2coeff(t.PAR);
	ids = 56:60; id = 58;

	[~,matF_gt{i},matM_gt{i},ind_gt{i}] = TD.estimate_3dtraj(im2double(t.Vk(:,:,:,ids)), gt_coeffs{i}(ids), params0.M, n, params_tbd3d);

	f1 = montage(matF_gt{i}(:,:,:,2*n+1:2*n+8),'Size',[1 8]); 
	f1 = f1.CData;
	m1 = montage(matM_gt{i}(:,:,:,2*n+1:2*n+8),'Size',[1 8]); m1 = m1.CData;
	imwrite(f1,['~/tmp/recon/tennis_F.png']);
	imwrite(m1,['~/tmp/recon/tennis_M.png']);

	r = size(f1,1)/2;
	fgt = [];
	for k = 1:8
		ind = (id-1)*n + k + 1;
		nma = sprintf('/mnt/lascar/rozumden/dataset/TbD_hs/tennis/%08d.png', ind);
		img = imread(nma);
		ps = gt_coeffs{i}{id}{k}(:,1) + gt_coeffs{i}{id}{k}(:,2);
		f = img(ps(2)-r:ps(2)+r,ps(1)-r:ps(1)+r,:);
		fgt = cat(2,fgt,f);
	end
	imwrite(fgt,['~/tmp/recon/tennis_FGT.png']);

	scl = 1;
	new_coeff = gt_coeffs{i}{id};
	for tmp = 1:n, new_coeff{tmp} = scl*new_coeff{tmp}; end
	% new_coeff = [gt_coeffs{i}{id}{1}(:,1) gt_coeffs{i}{id}{end}(:,1)];
	% new_coeff = {[new_coeff(:,1) new_coeff(:,2)-new_coeff(:,1)]*scl};
	inp = imresize(t.Vk(:,:,:,id),scl);
	Hall = myTrajRender(size2(inp), new_coeff, [0 1]);
	inp(logical(Hall>0)) = 0;
	inp(logical(cat(3,zeros(size(Hall)),Hall>0))) = 255;
	inp(logical(cat(3,zeros(size(Hall)),zeros(size(Hall)),Hall>0))) = 0;
	HM = conv2(Hall,ones(scl*size(matM_hs{i}(:,:,:,1))),'same');
	[yy,xx]=find(HM>0);

	% for tk = 1:n
	% 	% position = new_coeff{1}(:,1) + ((tk)/n)*new_coeff{1}(:,2);
	% 	position = new_coeff{tk}(:,1) + 1*new_coeff{tk}(:,2); 
	% 	position(2) = position(2)-10;
	% 	position(1) = position(1)-20;

	% 	inp = insertText(inp,position',tk,'BoxOpacity',0,'FontSize',22,'TextColor','white');
	% 	clr = [255 255 255];
	% 	% position = round(new_coeff{1}(:,1) + ((tk-0.5)/n)*new_coeff{1}(:,2));
	% 	position = new_coeff{tk}(:,1) + 0.5*new_coeff{tk}(:,2);
	% 	for ii = -1:1
	% 		for jj = -1:1
	% 			for ui = 1:3, inp(position(2)+ii,position(1)+jj,ui) = clr(ui); end
	% 		end
	% 	end
	% end
	inp = inp(min(yy(:)):max(yy(:)),min(xx(:)):max(xx(:)),:);
	imwrite(inp,['~/tmp/recon/tennis.png']);

	Hs = [];
	img = t.Vk(:,:,:,id);
	gap = 0.07;
	for k3 = 1:n
		H = evaluate_vcoeff(size2(t.Vk), gt_coeffs{i}{id}, [( ( (k3-1+gap) )/n) (k3/n)]);
		H = double(H);
		H = H ./ sum(H(:));
		Hs = cat(3, Hs, H);
	end
	Hs = Hs / sum(Hs(:));
	bgr = median(t.Vk,4);
	save(['~/tmp/tennis.mat'],'img','bgr','Hs','Hall','f1','m1','fgt');

end