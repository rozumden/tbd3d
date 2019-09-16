function [params, cfg] = get_params(fast_version)
if nargin < 1
	fast_version = false;
end
% Configurations of TbD (does not change performance)
cfg = [];
cfg.verbose = 1;
cfg.write = 0;
cfg.show = 1;
cfg.save_frames = 0;
cfg.method = 'TbD';
cfg.start_frame = 1;
cfg.end_frame = Inf;
cfg.get_visim = false; %% TODO: also should be changed in RobustDetector
cfg.extend_dataset = 0;

% Parameters of TbD (changes performance)
params = [];
params.BGR = [];
params.detec_type = 'fast';
params.fmod_max_deviation_from_sumh = [0.5 2.5];
params.resize = 1; %% 2/3
params.th = 5/255;
params.min_radius = 3;
params.expos = 0.8;
params.gm = 1; % gamma correction [0.48 0.47 0.48] or 1/2.2
params.crop_type = 'hmask'; % parallel or rotated or hmask
params.fitting_type = 'gradient_mask'; % gradient or dyn_prog or gradient_mask
params.em_cycles = 10;
params.tbd_alpha_h = 0.2; 
params.tbd_beta_h = 500; 
params.increase_edge = true;
params.do_tbd = false;
params.do_em = true;
params.linear_fit = false;
params.num_appear = 1;
params.s_th = 0.15;  
params.fmod_th = 0.4;
params.ext_factor = 2; 
params.reject_fit_on_edge = true;
params.do_stabilize = true;
params.check_on_M = false;
params.return_prediction = true;
params.use_template = 1;
params.alpha_F = 1; % F = (1-a)*F + a*F0
params.add_dp = 1;
params.tbd_on_fmod = 1;
params.distance_check = 1;
params.add_std = false;
params.long_term_template = 1;

params.loop_maxiter = 100;
params.loop_rel_tol = 5e-3;
params.cg_maxiter = 25;
params.seqransac_lin_max_samples = 7140;
params.ransac_quad_max_samples = 4000;

params.fast_version = fast_version;
params.low_contrast_threshold = 0.25;
params.apply_normalisation = true;

if fast_version
	params.loop_maxiter = 15;
	params.loop_rel_tol = 1e-2;
	params.cg_maxiter = 5;
	params.seqransac_lin_max_samples = 800;
	params.ransac_quad_max_samples = 1000;
	params.do_em = false;
end

