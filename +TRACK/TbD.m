function fr = TbD(im_orig, bgr_orig, fr0, params, cfg, template_pair)
fr = [];
coeff0 = cellfun(@(x)fliplr(x).', fr0.coeff, 'UniformOutput', false);
for c = 1:numel(coeff0)
	if ~isempty(fr0.bb)
		coeff0{c}(:,1) = coeff0{c}(:,1) + [fr0.bb(1:2)]' - 1;
	end
end
ext = params.ext_factor*fr0.Radius;

im_orig(im_orig < 0) = 0; im_orig(im_orig > 1) = 1;
bgr_orig(bgr_orig < 0) = 0; bgr_orig(bgr_orig > 1) = 1;
im = im_orig.^(1/params.gm);
bgr = bgr_orig.^(1/params.gm);

if isempty(fr0.Direction)
	[frmi, smi] = getFrame(im, bgr, fr0, coeff0, ext, -1, params, template_pair);
	[frpl, spl] = getFrame(im, bgr, fr0, coeff0, ext, 1, params, template_pair);
	fr = frpl;
	if smi > spl
		fr = frmi;
	end
else
	fr = getFrame(im, bgr, fr0, coeff0, ext, fr0.Direction, params, template_pair);
	if isempty(fr)
		fr = getFrame(im, bgr, fr0, coeff0, ext, -fr0.Direction, params, template_pair);
	end
	if isempty(fr) && params.return_prediction
		fr = getPredFrame(im, bgr, fr0, coeff0, ext, fr0.Direction, params, template_pair);
	end
end

if ~isempty(fr)
	%% get direction using closest points
	bb0 = fr0.bb; bb1 = fr.bb;
	if isempty(bb0), bb0 = fr0.mbb; end
	if isempty(bb1), bb1 = fr.mbb; end
	[fr.Direction, displacement] = FIT.get_dir(fr0.coeff, fr.coeff, bb0, bb1, fr0, fr);
	fr.Speed = fr.Length / fr.Radius;
	fr.SpeedFull = (fr.Length + displacement) / fr.Radius;

	fr.visim = [];
	if cfg.get_visim
		[~, T0] = FIT.trajPredictor(coeff0, fr.Direction, params.expos, [size(im,1) size(im,2)]);
		T0_c = T0(fr.bb(2):fr.bb(4),fr.bb(1):fr.bb(3),:);
		T0_c = T0_c / sum(T0_c(:));
		fr.visim = genvisim(T0_c, fr.T,fr.f,fr0.M,fr.h,fr.bgr_c,fr.im_c,zeros(size(fr.im_c)),2.2);
	end
end


function [fr, s, pts] = getFrame(im, bgr, fr0, coeff0, ext, direc, params, template_pair)
s_th = params.s_th;
fr = []; s = 0;
sz = [size(im,1) size(im,2)];

[bb, pts] = getPred(coeff0, ext, params.expos, sz, direc);
[~, ptsadd] = getPred(coeff0, ext, 1, sz, direc); pts = [pts; ptsadd];
if ~isempty(fr0.mbb)
	pts = TRACK.change_coor(fr0.mbb,pts(:,2),pts(:,1))';
end
mbb = [];
hmask = [];
if ~isempty(template_pair)
	template = template_pair.F;
	M = template_pair.M;
elseif params.use_template && strcmp(fr0.caseused, 'FMOd') %% it's FMO detector
	template = params.F;
	M = params.M;
else
	template = fr0.f;
	M = fr0.M;
end

if strcmp(params.crop_type, 'parallel')
	[im_c, bgr_c] = TRACK.make_crop_axis_parallel(im, bgr, bb);
	[~, ~, mbb_par] = TRACK.make_crop_axis(im, bgr, ext, pts);
elseif strcmp(params.crop_type, 'rotated')
	[im_c, bgr_c, mbb] = TRACK.make_crop_axis(im, bgr, ext, pts);
	bb = [];
	if ~isempty(fr0.bb) %% previous was not rotated
		wdir = mbb(:,2) - mbb(:,1);
		ang = radtodeg(atan(wdir(2)/wdir(1)));
		template = imrotate(fr0.f,ang,'bicubic','crop');
	end
	% [curve, curve_full, H] = FIT.fit_dp(h, true);
elseif strcmp(params.crop_type, 'hmask')
	T0 = zeros(size2(im));
	pts0 = round(pts); pts0(pts0 < 1) = 1;
	szf = [size(im,2) size(im,1)];
	szfm = repmat(szf, size(pts,1), 1);
	pts0(pts0 > szf) = szfm(pts0 > szf);
	T0(sub2ind(size(T0), pts0(:,2), pts0(:,1))) = 1;
	hmask = conv2(T0, double(diskMask([], round(ext))), 'same') > 0;
	[y,x] = find(hmask);
	bb = [min([x y]) max([x y])];
	[im_c, bgr_c] = TRACK.make_crop_axis_parallel(im, bgr, bb);
	hmask = hmask(bb(2):bb(4),bb(1):bb(3),:);
	T0_c = T0(bb(2):bb(4),bb(1):bb(3));
end
% if bbiou(bb, fr0.bb) > params.iou_nbr_frames_th, return; end

if min([size(im_c,1) size(im_c,2)]) < size(M,1), return; end
maxdiff = max(abs(im_c(repmat(T0_c>0,[1 1 3])) - bgr_c(repmat(T0_c>0,[1 1 3]))));
if maxdiff < 0.03
	return; 
end


if any(size2(im_c) < size2(template))
	return
end

[st, en] = getStEn(coeff0, direc, params, fr0, bb, mbb);
[h, f, m, T, coeff, len, s] = TRACK.TbD_loop(im_c, bgr_c, M, template, template, [], hmask, st, en, params, s_th, fr0.Speed, params.add_dp);
si = s;
dhmask = bwedge(hmask);

%% if small doesn't work, make it bigger
if TRACK.is_edge_bin(T > 0) > 0 || sum(sum(dhmask(T>0))) > 0 || s <= s_th(1) 
	hmask = conv2(T0, double(diskMask([], round(ext*3.5))), 'same') > 0;
	[y,x] = find(hmask);
	bb = [min([x y]) max([x y])];
	[im_c, bgr_c] = TRACK.make_crop_axis_parallel(im, bgr, bb);
	hmask = hmask(bb(2):bb(4),bb(1):bb(3),:);
	dhmask = bwedge(hmask);
	[st, en] = getStEn(coeff0, direc, params, fr0, bb, mbb);
	[h, f, m, T, coeff, len, s] = TRACK.TbD_loop(im_c, bgr_c, M, template, template, [], hmask, st, en, params, s_th, fr0.Speed, params.add_dp);
end

if params.add_dp && s > s_th(1)
	%% iterate with small mask
	for iterk = 1:3
		incr = 10;
		hmask = conv2(T, double(diskMask([], incr)), 'same') > 0;
		if sum(hmask(:)) == 0, s = 0; break; end
		if ~params.fast_version, h = []; end
		[h, f, m, T, coeff, len, s] = TRACK.TbD_loop(im_c, bgr_c, M, template, f, h, hmask, st, en, params, s_th, fr0.Speed, false);
		dhmask = bwedge(hmask);
		if s > 0 && params.reject_fit_on_edge && (TRACK.is_edge_bin(T > 0) > 0 || sum(sum(dhmask(T>0))) > 0)
			s = 0;
		end
		if s > 0, break; end
	end
end

if params.check_on_M && sum(m(:))/sum(fr0.M(:)) < 0.33
	s = 0;
end


if s > 0 && params.reject_fit_on_edge && (TRACK.is_edge_bin(T > 0) > 0 || sum(sum(dhmask(T>0))) > 0)
	s = 0;
end


len_diff = abs(fr0.Length - len)/fr0.Length; %% if too big (relative change)
len_small = len/fr0.Length; %% if too small
if fr0.Length < 10, len_small = 1; end

if s > s_th(1) && len_diff < 1 && len_small > 1/2
	fr = createFrame(im_c, bgr_c, h, f, T, coeff, bb, mbb, fr0, size(im));
	fr.M = (1-params.alpha_F)*double(m) + params.alpha_F*double(M);
	fr.f = (1-params.alpha_F)*f + params.alpha_F*template;
	fr.Length = len;
	fr.fittingScore = s;
	fr.hmask = hmask;
	fr.caseused='TbD';
	s = si;
end


function [fr] = getPredFrame(im, bgr, fr0, coeff0, ext, direc, params, template_pair)
sz = [size(im,1) size(im,2)];
[bb, pts] = getPred(coeff0, ext, params.expos, sz, direc);
p0 = pts(1,:);
if ~isempty(fr0.mbb)
	pts = TRACK.change_coor(fr0.mbb,pts(:,2),pts(:,1))';
end
if ~isempty(template_pair)
	template = template_pair.F;
	M = template_pair.M;	
elseif params.use_template && strcmp(fr0.caseused, 'FMOd') %% it's FMO detector
	template = params.F;
	M = params.M;
else
	template = fr0.f;
	M = fr0.M;
end
T0 = zeros(size2(im));
pts0 = round(pts); pts0(pts0 < 1) = 1;
szf = [size(im,2) size(im,1)];
szfm = repmat(szf, size(pts,1), 1);
pts0(pts0 > szf) = szfm(pts0 > szf);
T0(sub2ind(size(T0), pts0(:,2), pts0(:,1))) = 1;
hmask = conv2(T0, double(diskMask([], round(ext))), 'same') > 0;
[y,x] = find(hmask);
bb = [min([x y]) max([x y])];
[im_c, bgr_c] = TRACK.make_crop_axis_parallel(im, bgr, bb);
hmask = hmask(bb(2):bb(4),bb(1):bb(3),:);
T0_c = T0(bb(2):bb(4),bb(1):bb(3));
T0_c = conv2(T0_c, double(diskMask([], 2)), 'same') > 0;
T0_c = double(T0_c)/sum(double(T0_c(:)));

[coeff,T,~] = psffit(T0_c, params);
[~,~,len] = myTrajRender(size(T), cellfun(@(x)fliplr(x).', coeff, 'UniformOutput', false));
fr = createFrame(im_c, bgr_c, T0_c, template, T, coeff, bb, [], fr0, size(im));
fr.M = M;
fr.f = template;
fr.Length = len;
fr.fittingScore = 1;
fr.hmask = hmask;
fr.caseused = 'Pred';


function [st, en] = getStEn(coeff0, direc, params, fr0, bb, mbb)
pnts = FIT.trajPredictor(coeff0, direc, params.expos)';
if ~isempty(fr0.mbb)
	pnts = TRACK.change_coor(fr0.mbb,pnts(2,:)',pnts(1,:)');
	st = TRACK.change_coor_back(mbb, pnts(1,1), pnts(2,1))';
	en = TRACK.change_coor_back(mbb, pnts(1,end), pnts(2,end))';
else
	pnts = pnts - [bb(1:2)]'+1;
	st = pnts(:,1);
	en = pnts(:,end);
end

function [bb, pts] = getPred(coeff0, ext, expos, sz, direc)
pts = FIT.trajPredictor(coeff0, direc, expos);
p0 = round(min(pts) - ext);
p0(p0 < 1) = 1;
p1 = round(max(pts) + ext);
szf = fliplr(sz);
p1(p1 > szf) = szf(p1 > szf);
bb = [p0 p1];

function [bb] = getPredLocal(coeff0, ext, sz, direc, fr0)
if isempty(fr0.Start) || isempty(fr0.End)
	fr0.Start = coeff0{1}(:,1)';
	fr0.End = sum((coeff0{end})');
	fr0.Direction = 1;
end
pts = fr0.End;
if direc ~= fr0.Direction
	pts = fr0.Start;
end
ext1 = ext + fr0.Length;
p0 = round(pts - ext1);
p0(p0 < 1) = 1;
p1 = round(pts + ext1);
szf = fliplr(sz);
p1(p1 > szf) = szf(p1 > szf);
bb = [p0 p1];


function frm = createFrame(im_c, bgr_c, h, f, T, coeff, bb, mbb, fr0, sz)
frm = Frame();
frm.isFMO = true; 
frm.foundPrev = true;
frm.im_c = im_c; 
frm.bgr_c = bgr_c; 
frm.fmoCheck = true;
frm.f = f; 
frm.h = h; 
frm.T = T; 
frm.M = fr0.M; 
frm.coeff = coeff; 
frm.Radius = fr0.Radius;
frm.bb = bb;
frm.mbb = mbb;

TM = conv2(double(T),double(frm.M),'same');
[x,y] = find(T > 0);
[xx yy] = find(TM);

if isempty(mbb)
	frm.TrajectoryXY = [y'; x'] + bb(1:2)' - 1;
	frm.PixelList = [yy'; xx'] + bb(1:2)' - 1;
else
	frm.TrajectoryXY = TRACK.change_coor(mbb, x, y);
	frm.PixelList = TRACK.change_coor(mbb, xx, yy);
end

frm.PixelIdxList = sub2ind(sz,frm.PixelList(2,:),frm.PixelList(1,:));;
frm.TrajectoryXYF = frm.TrajectoryXY(2:-1:1,:);





