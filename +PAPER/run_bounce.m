[params, cfg] = EVAL.get_params(true);
[seq, folder] = EVAL.get_seq(0,'webcam');
params.th = 25/255;
params.use_template = false;
params.gm = 1;
seq(1).name = 'bounce_right';
for i = 1
	video = VID.VideoImg(folder, seq(i).name, 'webcam_img');

	% frame{i} = EVAL.tbd_main_loop(video, cfg, params);
	% VIS.show_all([], seq(i), frame{i});
end

n = 8;
gm = 1;
id = 75;
tm = 70:80;

bgr = [];
for tt = tm
	bgr = cat(4,bgr,video.get_frame(tt));
end
bgr = median(bgr,4);

img = video.get_frame(id);

p1 = [247.2367 236.8208]';
p2 = [220.9444 274.7468]';
p3 = [152.7708 215.8801]';
gt_cf = {[p1 p2-p1], [p2 p3-p2]};

inp = img;
Hall = myTrajRender(size2(inp), gt_cf, [0 1]);
inp(logical(cat(3,Hall>0))) = 0;
inp(logical(cat(3,zeros(size(Hall)),Hall>0))) = 255;
inp(logical(cat(3,zeros(size(Hall)),zeros(size(Hall)),Hall>0))) = 0;
HM = conv2(Hall,ones(size(M)),'same');
[yy,xx]=find(HM>0);
inp = inp(min(yy(:)):max(yy(:)),min(xx(:)):max(xx(:)),:);
imwrite(inp,['~/tmp/yellow_in.png']);

Hs = [];
gap = 0.07;
for k = 1:n
	H = evaluate_vcoeff(size2(img), gt_cf, [( ( (k-1+gap) )/n) (k/n)]);
	H = double(H);
	H = H ./ sum(H(:));
	Hs = cat(3, Hs, H);
end

maxiter =  10; 
lambda_R = 0; 
alpha_cross_f = 0; 
alpha_cross_m = 0; 
lambda_w = 1e-3; 
alpha_w = 2^-7;

r = 35;
M = double(diskMask([], r));
[Fs,Ms] = estimateFM_motion_pw(img, bgr, Hs, ones([size(M) 3 n]), repmat(M,[1 1 1 n]), [], [],...
	'alpha_f', alpha_w, 'alpha_m', 2^-8,  'lambda', lambda_w, 'maxiter', maxiter, ...
	'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6, ...
	'alpha_cross_f', alpha_cross_f, 'alpha_cross_m', alpha_cross_m, 'lambda_R', lambda_R);
Ms = permute(Ms,[1 2 4 3]);

f1 = montage(Fs,'Size',[1 size(Fs,4)]); f1 = f1.CData;
m1 = montage(Ms,'Size',[1 size(Fs,4)]); m1 = m1.CData;
fgt = f1;
save(['~/tmp/yellow_bounce.mat'],'img','bgr','Hs','Hall','f1','m1','fgt');
