function [tiou, tiou_mean] = gt_cost_iou_curves(curves, frms, PAR, r)
parts = 100;
% ef = 9/10; 
% ef = 11/12;
ef = 0.92;
tiou = zeros(1,numel(PAR));
for ki = 1:numel(curves)
	crv = curves(ki);
	if strcmp(crv.type,'connect')
		continue;
	end
	for ind = crv.fit_iv(1):crv.fit_iv(2)
		%% GT
		pars = PAR{ind};
		pars = pars(:,sum(pars) ~= 0);
		if size(pars,2) == 1, pars = [pars pars]; end
		coeff = {};
		for kk = 2:size(pars,2)
			st = pars(:,kk-1);
			en = pars(:,kk);
			vec = en-st;
			coeff = [coeff {fliplr([st'; vec'])}];
		end	
		[coeff, len, p2] = postprocc(coeff, [0; 0], parts);
		%%%%
		p1_raw = [];
		if ~isempty(frms{ind})
			[~, ~, p1_raw] = postprocc(frms{ind}.coeff, [frms{ind}.bb(1:2)]' - 1, parts);
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
		mn1 = mean(iou1);
		mn2 = mean(iou2);
		mn3 = mean(iou3);
		mn4 = mean(iou4);
		mn5 = mean(iou5);
		mn6 = mean(iou6);
		tiou(ind) = max([mn1 mn2 mn3 mn4 mn5 mn6]);
	end
end

tiou_mean = mean(tiou);

