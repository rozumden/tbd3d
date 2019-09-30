[params, cfg] = EVAL.get_params(true);
[seq, folder] = EVAL.get_seq(0,'webcam');
params.th = 25/255;
params.use_template = false;
params.gm = 1;
params.alpha_f = 0;
n = 8;

for i = 2
	video = VID.VideoImg(folder, seq(i).name, 'webcam_img');
	frame{i} = EVAL.tbd_main_loop(video, cfg, params);
	VIS.show_all([], seq(i), frame{i});
end


% [f_img,m_img] = VID.estimate_FM_mc(video, curves_all{i}, M, F);		
% [matF,matM] = TD.get_views(video,curves_all{i}, M, F);

[Fs,Ms] = TD.get_views_coeffs(frm);
[szs] = TD.estimate_3d(Fs,Ms,2);


% imgmean = zeros(size(video.get_frame(1)));
% for tt = tm, imgmean = imgmean + video.get_frame(tt); end
% imgmean = imgmean / numel(tm);


% [Fs,Ms,Hs] = VID.put_fmo_frame(fr,n,params.gm);
% [out] = VID.vis_rotation(fr.bgr_c, Hs, Fs, Ms, params.gm);
% [out] = VID.vis_rotation(fr.bgr_c, Hs(:,:,:,round(end/2)), Fs, Ms, params.gm);

% [szs, FF, MM] = EXPER.estimate_3d(frame{1}(tm), n);
% tms = linspace(0,numel(tm), numel(szs));
% szs2 = smooth(szs, 'rlowess');
% plot(tms, szs2);
