function [tiou, tiou_mean, len_gt] = gt_cost_iou_3d(frms0, PAR, method)
%% method: 1 - TIoU, 2 - TIoU-3D
if ~exist('method','var')
	method = 1;
end

frms = [frms0{:}];
frms = frms([frms.index] > 0);
parts = numel(PAR(1).R);

tiou = zeros(1,numel(frms0));
len_gt = zeros(1,numel(frms0));
usedind = [];
for ki = 1:numel(frms)
	frm = frms(ki);
	ind = frm.instance;
	usedind = [usedind ind];
	r = PAR(ind).R;

	pars = PAR(ind).POS;
	pars = pars(:,sum(pars) ~= 0);
	if isempty(pars) % FP
		costs_gt(ki) = 0;
	end
	if size(pars,2) == 1
		pars = [pars pars];
	end

	coeff = {};
	for kk = 2:size(pars,2)
		st = pars(:,kk-1);
		en = pars(:,kk);
		vec = en-st;
		coeff = [coeff {fliplr([st'; vec'])}];
	end	

	if ~isempty(frm.bb)
		[coeff1, len1, p1] = postprocc(frm.coeff, [frm.bb(1:2)]' - 1, parts);
	else
		[coeff1, len1, p1] = postprocc(frm.coeff, frm.mbb, parts);
	end
	[coeff, len, p2] = postprocc(coeff, [0; 0], parts);
	
	if method == 1
		iou1 = FIT.calciou(p1, p2, r);
		iou2 = FIT.calciou(fliplr(p1')', p2, r);
	elseif method == 2
		p1 = [p1 repmat(frm.Radius,parts,1)];
		p2 = [p2 r];
		iou1 = FIT.calciou3d(p1, p2, r);
		iou2 = FIT.calciou3d(fliplr(p1')', p2, r);
	else
		error('Method not defined');
	end

	mn1 = mean(iou1);
	mn2 = mean(iou2);
	tiou(frm.instance) = max(mn1,mn2);
	len_gt(frm.instance) = sum(len);
end

tiou_mean = mean(tiou);

