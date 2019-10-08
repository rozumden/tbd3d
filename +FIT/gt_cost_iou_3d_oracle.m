function [tiou, tiou_mean, gt_coeffs] = gt_cost_iou_3d_oracle(rads, PAR)
parts = numel(PAR(1).R);
tiou = zeros(1,numel(PAR));
gt_coeffs = {};
for ind = 1:numel(PAR)
	%% GT
	pars = PAR(ind).POS;
	r = PAR(ind).R;
	r = [r(:)];

	if numel(rads) == 1
		r_est = repmat(rads,parts,1);
	else
		if iscell(rads)
			r_est = rads{ind};
		else
			r_est = rads(ind).R;
		end
	end

	pars = pars(:,sum(pars) ~= 0);
	if size(pars,2) == 1, pars = [pars pars]; end
	coeff = {};
	for kk = 2:size(pars,2)
		st = pars(:,kk-1);
		en = pars(:,kk);
		vec = en-st;
		coeff = [coeff {fliplr([st'; vec'])}];
	end	
	[coeff, len, pnts] = postprocc(coeff, [0; 0], parts);
	gt_coeffs = [gt_coeffs {coeff}];

	p1 = [pnts r_est];
	p2 = [pnts r];
	iou1 = FIT.calciou3d(p1, p2, r);

	tiou(ind) = nanmean(iou1);
end

tiou_mean = mean(tiou);

