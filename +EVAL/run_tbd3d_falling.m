[params, cfg] = EVAL.get_params(true);

params_tbd3d = [];
params_tbd3d.do_hier = false;
params_tbd3d.f0_maxiter = 10;
params_tbd3d.maxiter = 5;
params_tbd3d.iter_smoothing = 4;
params_tbd3d.do_intervals = false;
params_tbd3d.lambda_R = 0;
% params_tbd3d.alpha_cross_f = 0;
% params_tbd3d.alpha_cross_m = 0;

[seq, folder] = EVAL.get_seq(0,'falling');
n = 10;
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

seq = seq(1:2);
Ms{2} = ones(50,25);
k = 27;
% gt_coeffs{2} = {{[182.8825  439.8902; -0.4013   47.3572]'}};
gt_coeffs{2} = {{[182.7362    0.8027; 456.5998   51.1030]}};
% [x,y] = getpts();
% gt_coeffs{2} = {{[x(1) y(1); x(2)-x(1) y(2)-y(1)]'}};
for i = 2:numel(seq)
	warning('off','all');
	disp(seq(i).name);
	params0 = params;
	seq(i).name = seq(i).name([1:end-11 end-8:end]);

	t = load(fullfile(folder, seq(i).name));
	Vk{i} = permute(t.V(:,end:-1:1,:,:),[2 1 3 4]);

	% params0.BGR = median(im2double(Vk{i}),4);
	% params0.M = double(diskMask([],r));
	% params0.F = repmat(params0.M,[1 1 3]);
	% video = VID.VideoMat(Vk{i});

	% Ms{i} = params0.M;

	if false
		timev = tic;
		frame{i} = EVAL.tbd_main_loop(video, cfg, params0);
		exectime(i) = toc(timev);
		fprintf('Sequence %s took %.3f sec.\n',seq(i).name, exectime(i));  

		VIS.show_all([], seq(i), frame{i});
	end

	% [expf, ~] = estimate_exposure(frame{i}); 
	% expf = 1;
	% frms = [frame{i}{:}];
	% r_est = max([frms.Radius]);
	% Ms{i} = double(diskMask([],r_est+5));

	% [curves{i}, frms] = FIT.sequence_fit(frame{i}, video.get_frame(3), expf, 6);

	% VIS.curve(curves{i});
	% VIS.curve(gt_coeffs{i});

	% [szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}), gt_coeffs{i}, Ms{i}, n, params_tbd3d);
	[szs_gt{i},matF_gt{i},matM_gt{i},ind{i}] = TD.estimate_3dtraj(im2double(Vk{i}(:,:,:,k:end)), gt_coeffs{i}, Ms{i}, n, params_tbd3d);

	% [szs_est{i},matF{i},matM{i},ind_est{i}] = TD.estimate_3dtraj(im2double(Vk{i}), curves{i}, Ms{i}, n, params_tbd3d);

	% [matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);

end

if false
	ki = k;
	save('~/projects/data/falling_res_10n_v2.mat','Ms','gt_coeffs','matF_gt','matM_gt','ki');
end

% EVAL.get_stats;
