function [tiou, tiou_mean, len_gt] = gt_cost_iou(frms0, PAR, r0)
if nargin < 3
	r0 = [];
end
frms = [frms0{:}];
frms = frms([frms.index] > 0);
parts = 100;

tiou = zeros(1,numel(frms0));
len_gt = zeros(1,numel(frms0));
usedind = [];
for ki = 1:numel(frms)
	frm = frms(ki);
	ind = frm.instance;
	usedind = [usedind ind];
	r = double(frm.Radius);
	if ~isempty(r0), r = r0; end

	pars = PAR{ind};
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
	
	iou1 = FIT.calciou(p1, p2, r);
	iou2 = FIT.calciou(fliplr(p1')', p2, r);

	mn1 = mean(iou1);
	mn2 = mean(iou2);
	tiou(frm.instance) = max(mn1,mn2);
	len_gt(frm.instance) = sum(len);
end

tiou_mean = mean(tiou);

