function [fr] = TbD_on_FMOd(frms, params, template_pair)
frs = [];

for frmi = 1:numel(frms)
	frm = frms(frmi);
	fr = [];
	if ~isempty(template_pair) %% reconstructed from prevous frames
		F = template_pair.F;
		M = template_pair.M;
	elseif params.use_template %% given from the beginning
		F = params.F;
		M = params.M;
	else %% whatever FMOd thinks
		F = frm.f;
		M = frm.M;
	end

	if any(size2(frm.im_c) < size2(F))
		return
	end

	im_c = frm.im_c.^(1/params.gm);
	bgr_c = frm.bgr_c.^(1/params.gm);

	% hmask = conv2(double(frm.T > 0), APP.gen_ball(round(1.5*frm.Radius)), 'same') > 0;
	params0 = params;
	% params0.loop_maxiter = min([100 2*params0.loop_maxiter]);

	% meanim = mean(mean(rgb2gray(im_c)));
	% if meanim < 0.1
	% 	params0.tbd_alpha_h = meanim;
	% 	params0.tbd_beta_h = 1000*params0.tbd_alpha_h;
	% end
	dif = sum(abs(im_c - bgr_c), 3); meandif = mean(dif(frm.Ti > 0));
	if meandif >= params.low_contrast_threshold
		params0.tbd_alpha_h = [];
		params0.tbd_beta_h = [];
	end
	
	[h, f, m, T, coeff, len, s] = TRACK.TbD_loop(im_c, bgr_c, M, F, F, [], [], [], [], params0, params.s_th, 0, 0, false);

	if sum(h(:)) < params.fmod_max_deviation_from_sumh(1) || sum(h(:)) > params.fmod_max_deviation_from_sumh(2)
		s = 0;
	end

	Ti = frm.Ti; di = bwdist(Ti>0); do = bwdist(T>0);
	dmean = (mean(di(T>0)) + mean(do(Ti>0)))/2;
	if dmean > frm.Radius + 5, s = 0; end

	if s >= params.fmod_th && TRACK.is_edge_bin(T > 0) == 0
		fr = Frame();
		fr.h = h; 
		fr.M = (1-params.alpha_F)*double(m) + params.alpha_F*double(M);
		fr.f = (1-params.alpha_F)*f + params.alpha_F*F;
		fr.im_c = frm.im_c;
		fr.bgr_c = frm.bgr_c;
		fr.coeff = coeff;
		fr.Length = len;
		fr.bb = frm.bb;
		fr.T = T;
		[x,y] = find(fr.T > 0);
		fr.TrajectoryXY = [y'; x'] + fr.bb(1:2)' - 1;
		fr.Radius = double(size(fr.M, 1)) / 2;
		fr.caseused='FMOd';
		fr.fittingScore = s;
		frs = [frs fr];
	end
end

fr = [];
if ~isempty(frs)
	[~,ind] = max([frs.fittingScore]);
	fr = frs(ind);
end

% if params.em_cycles == 213
% 	keyboard
% end