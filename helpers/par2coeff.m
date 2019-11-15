function [gt_coeffs] = par2coeff(PAR)
parts = 100;
gt_coeffs = {};
for ind = 1:numel(PAR)
	%% GT
	pars = PAR{ind};

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
				en = PAR{ind+1}(:,1);
			end
		else
			en = pars(:,kk+1);
		end
		vec = en-st;
		coeff = [coeff {fliplr([st'; vec'])}];
	end	
	[coeff, len, pnts] = postprocc(coeff, [0; 0], parts);
	gt_coeffs = [gt_coeffs {coeff}];
end
