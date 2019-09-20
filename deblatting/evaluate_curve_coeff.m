function [pnt] = evaluate_curve_coeff(curves, tind)
for ki = 1:numel(curves)
	crv = curves(ki);
	if strcmp(crv.type,'connect')
		continue;
	end
	if tind < crv.fit_iv(1) || tind > crv.fit_iv(2)+1,
		continue;
	end

	if strcmp(crv.type, 'bounce')
		%%% TO DO 
		pnt = evaluate_coeff(crv.coeff{1}, tind - crv.fit_iv(1));
	else
		pnt = evaluate_coeff(crv.coeff{1}, tind);
	end

	return;
end
