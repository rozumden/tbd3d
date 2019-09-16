function [seq, folder] = get_seq(n, dset)
if nargin < 2
	dset = 'TbD';
end
seq = [];

%% TODO: Change the folder for each datasets 

if strcmp(lower(dset), 'fmo')
	seq(numel(seq)+1).name = 'volleyball1.mp4';
	seq(end).resize = 1/3; seq(end).end_frame = 20;
	seq(numel(seq)+1).name = 'volleyball_passing.mp4';
	seq(end).resize = 2/3;
	seq(numel(seq)+1).name = 'darts1.mp4';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'darts_window1.mp4'; 
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'softball.avi'; 
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'william_tell.avi';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'tennis_serve_side.avi';
	seq(end).resize = 1; 
	seq(numel(seq)+1).name = 'tennis_serve_back.avi';
	seq(end).resize = 1; seq(end).end_frame = 75;
	seq(numel(seq)+1).name = 'tennis1.avi';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'hockey.avi';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'squash.avi';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'frisbee.mp4';
	seq(end).resize = 0.2;
	seq(numel(seq)+1).name = 'blue.mov';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'ping_pong_paint.mov';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'ping_pong_side.mp4';
	seq(end).resize = 1;
	seq(numel(seq)+1).name = 'ping_pong_top.mp4';
	seq(end).resize = 1;
	folder = '/mnt/lascar/rozumden/dataset';
elseif strcmp(lower(dset), 'tbd')
	% New sequence from GoPro HERO 7 and improved GT in VS_throw_floor in frames with shadow
	seq(numel(seq)+1).name = 'VS_badminton_white_GX010058-8.mat'; %% (equivalent to 30fps)
	seq(end).resize = 1; seq(end).start_ind = 46; seq(end).end_ind = 85; seq(end).short='badminton_white';
	seq(end).low_contrast = 1; seq(end).fps = 30; seq(end).non_uniform = 1;
	seq(numel(seq)+1).name = 'VS_badminton_yellow_GX010060-8.mat'; %% (equivalent to 30fps)
	seq(end).resize = 1; seq(end).start_ind = 69; seq(end).end_ind = 125; seq(end).short='badminton_yellow';
	seq(end).low_contrast = 1; seq(end).fps = 30; seq(end).non_uniform = 1;
	seq(end).linear_fit = true;
	seq(numel(seq)+1).name = 'VS_pingpong_GX010051-8.mat'; %% (equivalent to 30fps)
	seq(end).resize = 1; seq(end).start_ind = 28; seq(end).end_ind = 85; seq(end).short='pingpong';
	seq(end).low_contrast = 1; seq(end).fps = 30; seq(end).gt_offset = 14;
	seq(end).linear_fit = true;
	seq(numel(seq)+1).name = 'VS_tennis_GX010073-8.mat'; %% (equivalent to 30fps)
	seq(end).resize = 1; seq(end).start_ind = 31; seq(end).end_ind = 68; seq(end).short='tennis';
	seq(end).do_em = false; seq(end).linear_fit = true;
	seq(end).low_contrast = 1; seq(end).fps = 30;
	seq(numel(seq)+1).name = 'VS_volleyball_GX010068-12.mat'; %% (equivalent to 20fps)
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 41; seq(end).short='volleyball';
	seq(end).low_contrast = 0; seq(end).fps = 20; seq(end).gt_offset = 17;  seq(end).gt_offset_hs = -5;
	seq(end).non_uniform = 1;
	seq(numel(seq)+1).name = 'VS_throw_floor-gc-8_newGT.mat'; %% (better GT in areas with shadows but all maybe slightly shifted upwards)
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 40; seq(end).short='throw_floor';
	seq(end).low_contrast = 0; seq(end).fps = 30;
	seq(numel(seq)+1).name = 'VS_throw_soft-gc-4_newGT.mat'; %% (better GT)
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 60; seq(end).short='throw_soft';
	seq(end).low_contrast = 1; seq(end).fps = 60;
	seq(end).non_uniform = 1;
	seq(numel(seq)+1).name = 'VS_throw_tennis-gc-8.mat';
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 45; seq(end).short='throw_tennis';
	seq(end).low_contrast = 0; seq(end).fps = 30;
	seq(numel(seq)+1).name = 'VS_roll_golf-gc-12.mat';
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = []; seq(end).short='roll_golf';
	seq(end).pad_pcent = 0.5;  seq(end).fps = 20;
	seq(end).low_contrast = 0; seq(end).non_uniform = 1;
	seq(end).do_em = false; seq(end).linear_fit = true;
	seq(numel(seq)+1).name = 'VS_fall_cube-gc-8.mat';
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 20; seq(end).short='fall_cube';
	seq(end).low_contrast = 1; seq(end).fps = 30; seq(end).non_uniform = 1;
	seq(numel(seq)+1).name = 'VS_hit_tennis-gc-8.mat';
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = 30; seq(end).short='hit_tennis';
	seq(end).low_contrast = 0; seq(end).fps = 30; seq(end).pad_pcent = 0.2;
	seq(numel(seq)+1).name = 'VS_hit_tennis2-gc-8.mat';
	seq(end).resize = 2/3; seq(end).start_ind = 1; seq(end).end_ind = []; seq(end).short='hit_tennis2';
	seq(end).low_contrast = 0; seq(end).fps = 30; seq(end).pad_pcent = 0.2;
	folder = '/mnt/lascar/rozumden/dataset/TbD_GC';
elseif strcmp(lower(dset), 'atp')
	% 1
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 1; seq(end).end_frame = 67; % 67
	seq(end).x1 = 200; seq(end).x2 = 700; 
	seq(end).hit = [1720 375];

	% 2
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 68; seq(end).end_frame = 137; % 137
	seq(end).x1 = 400; seq(end).x2 = 900; 
	seq(end).hit = [1680 500];

	% 3
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 138; seq(end).end_frame = 199;
	seq(end).hit = [1680 500];

	% 4
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 200; seq(end).end_frame = 274;
	seq(end).hit = [1680 500];

	% 5
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 275; seq(end).end_frame = 356;

	% 6
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 357; seq(end).end_frame = 386;

	% 7
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 387; seq(end).end_frame = 420;

	% 8
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 421; seq(end).end_frame = 498;

	% 9
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 499; seq(end).end_frame = 565;

	% 10
	seq(numel(seq)+1).name = 'atp_serves.avi';
	seq(end).start_frame = 566; seq(end).end_frame = 655;
	folder = '/mnt/lascar/rozumden/dataset';
else
	error('Dataset not found');
end	

if nargin > 0 && n > 0
	seq = seq(n);
end
