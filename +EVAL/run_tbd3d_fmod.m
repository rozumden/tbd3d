
[seq, folder] = EVAL.get_seq(0,'TbD-3D');
resz = 0.5;
verbose = false;

% Vk = repmat({[]}, 1, numel(seq));

for n = [8 4 2 1]

	% for i = 5
	% for i = 1:numel(seq)
	if isempty(gcp('nocreate')), parpool(numel(seq)); end
	parfor i = 1:numel(seq)
		warning('off','all');
		disp(seq(i).name);

		% if isempty(Vk{i}) 
			t = load(fullfile(folder, seq(i).name));
			if ~isempty(seq(i).crop), t.V = t.V(1:seq(i).crop(1),1:seq(i).crop(2),:,:); end
			[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
			Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
			Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);
		% end

		BGR = median(im2double(Vk{i}),4);
		video = VID.VideoMat(Vk{i});
		
		clf; image(video.get_frame(1));

		timev = tic;
		frame{i} = cell(1,video.length());
		detec = DET.SimpleDetector('BGR',BGR,'noise',35/255,'linear_fit',false,'linear_fit_init',false,'do_deblur',false);
		% try
			for k = 1:video.length()
				timev2 = tic;
				frame{i}{k} = detec.detect(video.get_frame(k));
				if k > 1
					frame{i}{k}.Direction = FIT.get_dir(frame{i}{k-1}.coeff, frame{i}{k}.coeff, frame{i}{k-1}.bb, frame{i}{k}.bb, frame{i}{k-1}, frame{i}{k});
				end
				frame{i}{k}.show; drawnow;
				if verbose, fprintf('Frame %d/%d, %.3f sec \n',k,video.length(),toc(timev2)); end
			end

			exectime(i) = toc(timev);
			if verbose, fprintf('Sequence %s took %.3f sec.\n',seq(i).name, exectime(i)); end

			% [expf, ~] = estimate_exposure(frame{i}); 
			expf = 1;
			frms = [frame{i}{:}];
			r_est = max([frms.Radius]);
			Ms{i} = double(diskMask([],r_est+5));
			
			[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});

			VIS.show_all([], seq(i), frame{i});
			VIS.curve(gt_coeffs{i});
			
			[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
			fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

			timev = tic;
			[curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6, true, 3);
			exectime_nc(i) = toc(timev);
			if verbose, fprintf('[NC] Sequence %s took %.3f sec.\n',seq(i).name, exectime_nc(i)); end

			[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
			fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));

			VIS.curve(curves{i});
		% catch err
		% 	disp(['Error: ' err]);
		% end
		[~,seqname,~] = fileparts(seq(i).name);
		C = strsplit(seqname,'_');
		short = C{end};
		saveas(gcf,['~/tmp/' short '_res' int2str(n) '.png']);
	end

	EVAL.get_stats;

	save(['~/projects/data/TbD-3D-n' int2str(n) '.mat'],'frame','curves','gt_coeffs');
end