function [tiou, tiou_mean] = gt_cost_iou_3d_curves(curves, frms, PAR, method, r3d, inds)
%% method: 1 - TIoU, 2 - TIoU-3D
if ~exist('method','var')
	method = 1;
end

if ~exist('inds','var')
	inds = [];
end
if ~isempty(inds)
	rr = r3d;
	if iscell(r3d)
		rr = [r3d{:}]; 
	end
	rr = [rr(:)];
end

parts = numel(PAR(1).R);
% ef = 9/10; 
% ef = 11/12;
% ef = 0.92;
ef = 1;

frms_all = [frms{:}];
r_est = frms_all(1).Radius;

tiou = zeros(1,numel(PAR));
for ki = 1:numel(curves)
	crv = curves(ki);
	if strcmp(crv.type,'connect')
		continue;
	end
	for ind = crv.fit_iv(1):crv.fit_iv(2)
		%% GT
		pars = PAR(ind).POS;
		r = PAR(ind).R;
		r = [r(:)];
		
		pars = pars(:,sum(pars) ~= 0);
		if size(pars,2) == 1, pars = [pars pars]; end
		coeff = {};
		vec = [];
		for kk = 1:size(pars,2)
			st = pars(:,kk);
			if kk == size(pars,2)
				if ind == numel(PAR)
					en = st + vec;
				else
					en = PAR(ind+1).POS(:,1);
				end
			else
				en = pars(:,kk+1);
			end
			vec = en-st;
			coeff = [coeff {fliplr([st'; vec'])}];
		end	
		[coeff, len, p2] = postprocc(coeff, [0; 0], parts);
		%%%%
		p1_raw = [];
		if ~isempty(frms{ind})
			[~, ~, p1_raw] = postprocc(frms{ind}.coeff, [frms{ind}.bb(1:2)]' - 1, parts);
			r_est = frms{ind}.Radius;
		end
		if strcmp(crv.type,'joint') || strcmp(crv.type,'prediction')
			[~, ~, p1] = postprocc(crv.coeff, [], parts, [ind ind+ef], false);
			[~, ~, p1_ref] = postprocc(crv.coeff, [], parts, [ind+(1-ef) ind+1], false);
		elseif strcmp(crv.type,'bounce')
			ilen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
			one_part_len = round(parts/ef);
			if ilen == 1
				[~, ~, p1] = postprocc(crv.coeff, [], parts, [0 1], false);
				[~, ~, p1_ref] = postprocc(crv.coeff, [], one_part_len, [0 1], false);
				p1_ref = p1_ref(1:parts,:);
			else
				parts0 = one_part_len * ilen;
				[~, ~, p1_ext] = postprocc(crv.coeff, [], parts0, [0 1], false);
				crvind = ind - crv.fit_iv(1);
				st_ind = one_part_len*crvind + 1;
				p1 = p1_ext(st_ind:st_ind+parts-1,:);
				p1_ref = p1_ext(st_ind+one_part_len-parts:st_ind+one_part_len-1,:);
			end
		end
		
		if method == 1
			iou1 = FIT.calciou(p1, p2, r);
			iou2 = FIT.calciou(fliplr(p1')', p2, r);
			iou3 = FIT.calciou(p1_ref, p2, r);
			iou4 = FIT.calciou(fliplr(p1_ref')', p2, r);
			if ~isempty(p1_raw)
				iou5 = FIT.calciou(p1_raw, p2, r);
				iou6 = FIT.calciou(fliplr(p1_raw')', p2, r);
			else
				iou5 = 0; iou6 = 0;
			end
		elseif method == 2
			if exist('r3d','var')
				if isempty(inds)
					r_use = r3d{ind};
				else
					xq = linspace(ind, ind+1, parts+1); xq = xq(1:end-1);
					r_use = interp1(inds,rr,xq)';
				end
			else
				r_use = repmat(r_est,parts,1);
			end
			p1 = [p1 r_use];
			p1_ref = [p1_ref r_use];
			p1_raw = [p1_raw r_use];
			p2 = [p2 r];
			iou1 = FIT.calciou3d(p1, p2, r);
			iou2 = FIT.calciou3d(fliplr(p1')', p2, r);
			iou3 = FIT.calciou3d(p1_ref, p2, r);
			iou4 = FIT.calciou3d(fliplr(p1_ref')', p2, r);
			if ~isempty(p1_raw)
				iou5 = FIT.calciou3d(p1_raw, p2, r);
				iou6 = FIT.calciou3d(fliplr(p1_raw')', p2, r);
			else
				iou5 = 0; iou6 = 0;
			end
		else
			error('Method not defined');
		end

		mn1 = nanmean(iou1);
		mn2 = nanmean(iou2);
		mn3 = nanmean(iou3);
		mn4 = nanmean(iou4);
		mn5 = nanmean(iou5);
		mn6 = nanmean(iou6);
		tiou(ind) = max([mn1 mn2 mn3 mn4 mn5 mn6]);
	end
end

tiou_mean = mean(tiou);

