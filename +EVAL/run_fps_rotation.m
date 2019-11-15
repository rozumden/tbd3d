resz = 0.5;
do_save = false;
do_parallel = true;

[params, cfg] = EVAL.get_params(true);
[seq, folder] = EVAL.get_seq(0,'TbD-3D');


[seq, folder] = EVAL.get_seq(0,'TbD-3D');

fc = [1.9 1 1.8];
WB = [2 1 2]; gamma_coef = 0.4;

rot_averages = struct('err_or_u',[],'err_u',[],'err_or_rot',[],'err_rot',[],'alla',[]);

if do_parallel
	seq_id = [1:numel(seq)];
	ns = [8 4 2 1];
else
	seq_id = 7;
	ns = [1]; 
end

for n = ns
	matF_hs_cell = repmat({[]}, 1, numel(seq));
	matF_gt_cell = repmat({[]}, 1, numel(seq));
	matF_cell = repmat({[]}, 1, numel(seq));

	load(['~/projects/data/TbD-3D-n' int2str(n) '.mat']);
	load(['~/projects/data/TbD-3D-n' int2str(n) '_post.mat']);
	
	if do_parallel && isempty(gcp('nocreate')), parpool(numel(seq)); end
	fname = ['~/projects/data/TbD-3D-n' int2str(n) '_matF.mat'];
	if exist(fname,'file')
		load(fname);
	else
		parfor i = seq_id
		% for i = seq_id
			disp(seq(i).name);
			warning('off','all');
		
			t = load(fullfile(folder, seq(i).name));
			[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
			Vk{i}(:,:,1,:) = fc(1).*Vk{i}(:,:,1,:);
			Vk{i}(:,:,3,:) = fc(3).*Vk{i}(:,:,3,:);
			[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);
			
			rr = [PAR{i}.R];
			ind_r = linspace(1,numel(frame{i})+1,numel(rr));
			matF_hs_cell{i} = matF_rescale(matF_hs{i}, [rr(:)]);
			% matF_gt_cell{i} = matF_rescale(matF_gt{i}, szs_gt{i});
			% matF_cell{i} = matF_rescale(matF{i}, szs_est{i});
			rr_gt = interp1(ind_r, [rr(:)], ind_gt{i},'linear', 'extrap');
			rr_est = interp1(ind_r, [rr(:)], ind_est{i},'linear', 'extrap');
			matF_gt_cell{i} = matF_rescale(matF_gt{i}, rr_gt);
			matF_cell{i} = matF_rescale(matF{i}, rr_est);
		end
		if do_save
			save(fname,'matF_gt_cell','matF_cell','matF_hs_cell','ind_gt','ind_est');
		end
	end

	EVAL.run_tbd3d_rotation

	rot_averages.err_or_u = [rot_averages.err_or_u mean(err_or_u)];
	rot_averages.err_u = [rot_averages.err_u mean(err_u)];
	rot_averages.err_or_rot = [rot_averages.err_or_rot mean(err_or_rot)];
	rot_averages.err_rot = [rot_averages.err_rot mean(err_rot)];
	aalstr.err_or_u = err_or_u;
	aalstr.err_u = err_u;
	aalstr.err_or_rot = err_or_rot;
	aalstr.err_rot = err_rot;
	rot_averages.alla = [rot_averages.alla aalstr];
end

if do_save
	save('~/projects/data/rotation_averages.mat','rot_averages');
	% save('~/projects/data/TbD-3D-rotations-sliding.mat','rot_averages');
	% save('~/projects/data/TbD-3D-rotations.mat', 'err_or_u','err_u','err_or_rot','err_rot');
end

