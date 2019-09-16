function [h, f, m, T, coeff, len, s] = TbD_loop(im_c, bgr_c, M, template, f, h, hmask, st, en, params, s_th, speed, add_dp, do_normalize)
if ~exist('do_normalize','var') 
	do_normalize = true;
end
if ~exist('speed','var') 
	speed = [];
end
if ~exist('add_dp','var') 
	add_dp = false;
end
coeff = []; m = M; T = h;
len = 0;
s = 0;

if ~isempty(hmask)
	dif = sum(abs(im_c - bgr_c), 3); 
	tr0 = bwmorph(hmask,'thin',inf);
	meandif = mean(dif(tr0 > 0));
	if meandif >= params.low_contrast_threshold
		params.tbd_alpha_h = [];
		params.tbd_beta_h = [];
	end
end

if ~isempty(hmask)
	roi = boundingBox(hmask);
	fsize = size(f);
	pad = ceil((fsize(1:2)-1)/2);
	roi = roi + [-pad(1), +pad(1), -pad(2), +pad(2)];
	hsize_orig = size2(im_c);
	roi = min(max(roi, [1 1 1 1]), hsize_orig([1 1 2 2])); 
	roisz = roi([2 4]) - roi([1 3]) + 1;
	if any(roisz < size(f,1))
		return
	end
end

if strcmp(params.fitting_type, 'gradient_mask')
	maxiter = params.loop_maxiter;
	rel_tol = params.loop_rel_tol;

	if isfield(params, 'tbd_alpha_h') && isfield(params, 'tbd_beta_h') && ...
		~isempty(params.tbd_alpha_h) && ~isempty(params.tbd_beta_h)
		[h,f,m,~] = blindLoop_FMH_TbD(im_c, bgr_c, f, M, h, hmask, template, 'maxiter', maxiter, ...
		    'rel_tol', rel_tol, 'alpha_h', params.tbd_alpha_h,'beta_h', params.tbd_beta_h, ...
		    'cg_maxiter', params.cg_maxiter);
	else
		[h,f,m,~] = blindLoop_FMH_TbD(im_c, bgr_c, f, M, h, hmask, template, 'maxiter', maxiter,...
			'rel_tol', rel_tol, 'cg_maxiter', params.cg_maxiter);
	end
	if add_dp
		htemp = h; htemp(~hmask) = -Inf;
		[curve, hh] = FIT.fit_dp(htemp, st, en, false);
		T = zeros(size(h)); 
		T(sub2ind(size(h),curve(2,:), curve(1,:))) = 1; %% make new psf from DP fit 
		T = T / sum(T(:)); T = imgaussfilt(T,1); 
		s = sum(h(T>0));
		coeff = [];
	else
		[coeff,T,val] = psffit(h, params);
		if ~isempty(coeff{1})
			[~,~,len] = myTrajRender(size(h), cellfun(@(x)fliplr(x).', coeff, 'UniformOutput', false));
		else
			coeff{1} = {[1 1]};
		end
		T0s = conv2(double(T>0),ones(3),'same')>0;
		s = sum(h(T0s));
		if val < 0.4
			s = 0;
		end
	end
elseif strcmp(params.fitting_type, 'gradient') || (~isempty(speed) && speed < 1.0)
	[h,f,~,hi] = blindLoop_FH_TbD(im_c, bgr_c, M, f, h, hmask);
	[T, coeff, s, len] = FIT.fitting(h, s_th(1), 1, params);
elseif strcmp(params.fitting_type, 'dyn_prog')
	[h,f,~,hi] = blindLoop_FH_TbD(im_c, bgr_c, M, f, h, hmask);
	htemp = h; htemp(~hmask) = -Inf;
	[curve, hh] = FIT.fit_dp(htemp, st, en, false);
	T_pre = zeros(size(h)); 
	T_pre(sub2ind(size(h),curve(2,:), curve(1,:))) = 1; %% make new psf from DP fit 
	T_pre = T_pre / sum(T_pre(:)); T_pre = imgaussfilt(T_pre,1); %% suitable for Honza fitting
	[T, coeff, s, len] = FIT.fitting(T_pre, s_th(1), 1, params); %% fiting on DP
	T0s = conv2(double(T>0),ones(3),'same')>0;
	s = sum(h(T0s));
	% [~, len, p1] = postprocc(coeff_grad, frm.mbb, parts);
else
	error('Fitting type unknown');
	% hext = conv2(T_pre, ones(3), 'same');
	% [T, coeff, s, len] = FIT.getpoint(h, hext, s_th(1), max(hext(:)));
end
T = T/sum(T(:));
if do_normalize
	h = h/sum(h(:));
end
