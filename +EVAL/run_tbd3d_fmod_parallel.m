
[seq, folder] = EVAL.get_seq(0,'TbD-3D');
n = 8;
resz = 1;

for i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name);

	if isempty(Vk{i})
		t = load(fullfile(folder, seq(i).name));
		if ~isempty(seq(i).crop), t.V = t.V(1:seq(i).crop(1),1:seq(i).crop(2),:,:); end
		[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
		Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
		Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);
	end
	if isempty(BGR{i})
		BGR{i} = median(im2double(Vk{i}),4);
	end
	frame{i} = cell(1,size(Vk{i},4));
end

maxfrms = max(cellfun(@(x) size(x,4), Vk));
if isempty(gcp('nocreate')), parpool(maxfrms); end

framepar = cell(1,maxfrms);
parfor k = 1:maxfrms
	framepar{k} = cell(1,numel(seq));
	for i = 1:numel(seq)
		detec = DET.SimpleDetector('BGR',BGR{i},'noise',20/255,'linear_fit',true);
		if k > size(Vk{i},4), continue; end
		framepar{k}{i} = detec.detect(im2double(Vk{i}(:,:,:,k)));
	end
	disp(['Worker ' int2str(k) ' finished']);
end

for i = 1:numel(seq)
	for k = 1:size(Vk{i},4)
		frame{i}{k} = framepar{k}{i};
		if k > 1
			frame{i}{k}.Direction = FIT.get_dir(frame{i}{k-1}.coeff, frame{i}{k}.coeff, frame{i}{k-1}.bb, frame{i}{k}.bb, frame{i}{k-1}, frame{i}{k});
		end
	end
	VIS.show_all([], seq(i), frame{i});
end

parfor i = 1:numel(seq)
	warning('off','all');
	expf = 1;
	
	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	timev = tic;
	[curves{i}, frms] = FIT.sequence_fit(frame{i}, Vk{i}(:,:,:,3), expf, 6);
	exectime_nc(i) = toc(timev);
	fprintf('[NC] Sequence %s took %.3f sec.\n',seq(i).name, exectime_nc(i));  

	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));
end

EVAL.get_stats;
