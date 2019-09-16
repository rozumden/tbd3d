function [frame,frameall] = single_tbd_only(video, cfg, params)

sz = [video.height() video.width()];
vsz = video.length();
resz = params.resize;
bgr = [];
if ~isempty(params.BGR)
	bgr = imresize(params.BGR, resz);
end


detec = DET.FastDetector(sz*resz, 'BGR', bgr, 'noise', params.th, 'min_radius', params.min_radius, 'do_stabilize', params.do_stabilize);


if cfg.show, figure(1); clf; end
time = [];
if ~isempty(cfg.end_frame) && cfg.end_frame < vsz
	vsz = cfg.end_frame;
end
tracker = [];
IM0 = [];
template_pair = [];

for k = cfg.start_frame:vsz
	kind = k-params.num_appear+1;
	kind = max([cfg.start_frame kind]);
	try
		IM = imresize(video.get_frame(k), resz); % next step 
		IMproc = imresize(video.get_frame(kind), resz);
	catch
		disp('Error while reading! Aborting...');
		break;
	end
	tic; 

	frame{kind} = DET.init_det(detec, IM);
	used = 'FMOd';
	
	if isempty(frame{kind}), used = '-'; end
	
	frmtr = [];
	frmprev = [];
	if kind > cfg.start_frame && ~isempty(frame{kind-1}) % TbD
		background = detec.B;
		frmprev = frame{kind-1};

		if ~isempty(frmprev.f)
			frmtr = TRACK.TbD(IMproc, background, frmprev, params, cfg, template_pair);
		end
	end

	%% if there was some FMOd and just Pred
	if params.tbd_on_fmod && ~isempty(frame{kind}) && (isempty(frmtr) || strcmp(frmtr.caseused,'Pred'))
		newfrm = TRACK.TbD_on_FMOd(frame{kind}, params, template_pair);
		if ~isempty(newfrm)
			used = 'FMOd->TbD';
		else
			used = '- (TbD)';
		end
		frame{kind} = newfrm;
	end

	%% If it was just a prediction (thus failed) and no FMOd and low speed
	if params.add_std && ~isempty(frmprev) && isempty(frame{kind}) && (isempty(frmtr) || (strcmp(frmtr.caseused,'Pred') && frmtr.Speed <= 1))
		[frmtr, tracker] = TRACK.sota(frmprev, uint8(255*IM0), uint8(255*IMproc), tracker, params.std_version);
		frmtr.caseused = 'Std';
	end
	
	use_tracker = true;
	if ~isempty(frmtr) && ~strcmp(frmtr.caseused,'TbD') && ~isempty(frame{kind})
		use_tracker = false;
	end

	if ~isempty(frmtr) && use_tracker
		used = frmtr.caseused;
		if ~strcmp(frmtr.caseused, 'TbD')
			used = ['TbD->' frmtr.caseused];
		end			
		frame{kind} = frmtr;
	end
	%%%%

	if params.distance_check && kind > cfg.start_frame && ~isempty(frame{kind}) 
		indfr0 = find(cellfun(@(x) ~isempty(x) && ~strcmp(x.caseused, 'FMOd') && ~strcmp(x.caseused, 'Pred'), frame(1:kind-1)));
		if ~isempty(indfr0)
			indfr0 = indfr0(end);
			fr0 = frame{indfr0};
			d = norm(fr0.bb(1:2) - frame{kind}.bb(1:2));
			dind = kind - indfr0;
			th = max([3*fr0.Length+6*fr0.Radius 10*fr0.Radius]);
			if d > dind*th
				if isempty(frmtr)
					frame{kind} = [];
					used = '- (too far)';
				else
					frame{kind} = frmtr;
					used = ['TbD->' frmtr.caseused];
				end
			end
		end
	end

	if params.fast_version && ~isempty(frame{kind}) && strcmp(frame{kind}.caseused,'Pred') && ...
		~isempty(frame{kind-1}) && strcmp(frame{kind-1}.caseused,'Pred') && ...
		~isempty(frame{kind-2}) && strcmp(frame{kind-2}.caseused,'Pred') 
		frame{kind} = [];
		used = ['- (pred)'];
	end

	if ~isempty(frame{kind}) && strcmp(frame{kind}.caseused, 'TbD') && params.long_term_template
		template_pair.F = frame{kind}.f;
		template_pair.M = frame{kind}.M;
	end

	if ~isempty(frame{kind}) && ~strcmp(frame{kind}.caseused, 'Std')
		tracker = [];
	end		
	IM0 = IMproc;

	time = [time toc];

	if cfg.show
		VIS.show_main(IMproc,frame,kind); 
	end
	if cfg.verbose, fprintf('Frame %d/%d, %.3f sec, dets %d: %s.\n',kind,vsz,time(end),numel(frame{kind}),used); end 
end

if cfg.verbose, fprintf('Mean fps %.1f\n', 1/mean(time)); end

