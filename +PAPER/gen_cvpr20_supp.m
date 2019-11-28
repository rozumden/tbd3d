[params, cfg] = EVAL.get_params(true);
[seq, folder] = EVAL.get_seq(0,'TbD-3D');
i = 7;
frmi = 1;

bgr = median(im2double(Vk{i}),4);
img = im2double(Vk{i}(:,:,:,frmi));

Hscell = {};
gap = 0.0;
for ni = [0:3]
	n = 2^ni;
	Hs = [];
	sten = [];
	for k = 1:n
		H = evaluate_vcoeff(size2(img), gt_coeffs{i}{frmi}, [( ( (k-1+gap) )/n) (k/n)]);
		H = double(H);
		H = H ./ sum(H(:));
		Hs = cat(3, Hs, H);
	end
	Hs = Hs / sum(Hs(:));
	Hscell{ni+1} = Hs;
end

starts = {};
ends = {};
for k = 1:8
	starts{k} = gt_coeffs{i}{frmi}{k}(:,1);
	ends{k} = gt_coeffs{i}{frmi}{k}(:,1) + gt_coeffs{i}{frmi}{k}(:,2);
end


maxiter =  10; 
lambda_R = 1e-2; 
alpha_cross_f = 2^-12; 
alpha_cross_m = 2^-12; 
lambda_w = 1e-3; 
alpha_w = 2^-12;

M = Ms{i};

Fscell = {};
Mscell = {};
for ni = [0:3]
	n = 2^ni;
	Ft = ones([size(M) 3 n]);
	Mt = repmat(M,[1 1 1 n]);
	if n > 1
		Ft = Fs(:,:,:,repelem(1:size(Fs,4),2));
		Mt = Ms(:,:,repelem(1:size(Fs,4),2));
	end

	[Fs,Ms] = estimateFM_motion_pw(img, bgr, Hscell{ni+1}, Ft, Mt, [], [],...
		'alpha_f', alpha_w, 'alpha_m', 2^-12,  'lambda', lambda_w, 'maxiter', maxiter, ...
		'rel_tol', 0, 'cg_maxiter', 50, 'cg_tol', 1e-6, ...
		'alpha_cross_f', alpha_cross_f, 'alpha_cross_m', alpha_cross_m, 'lambda_R', lambda_R);
	Ms = permute(Ms,[1 2 4 3]);

	% f1 = montage(Fs,'Size',[1 size(Fs,4)]); f1 = f1.CData;
	% m1 = montage(Ms,'Size',[1 size(Fs,4)]); m1 = m1.CData;
	Fscell{ni+1} = Fs;
	Mscell{ni+1} = Ms;
end

save(['~/volleyball_hierarchy.mat'],'img','bgr','Hscell','Fscell','Mscell','starts','ends');
