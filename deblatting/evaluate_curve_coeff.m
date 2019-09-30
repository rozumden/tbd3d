function [pnt] = evaluate_curve_coeff(curves, tind)
parts = 100;
for ki = 1:numel(curves)
	crv = curves(ki);
	if strcmp(crv.type,'connect')
		continue;
	end
	if tind < crv.fit_iv(1) || tind > crv.fit_iv(2)+1,
		continue;
	end

	if strcmp(crv.type, 'bounce')
		if numel(crv.coeff) ~= 2, error('Should be a polynom!'); end
		[~, len, ~] = postprocc(crv.coeff, [], parts, [0 1], false);
		crvlen = crv.fit_iv(2) - crv.fit_iv(1) + 1;
		brkpnt = (len / sum(len)) * crvlen;
		tind0 = tind - crv.fit_iv(1);
		if tind0 <= brkpnt(1)
			pnt = evaluate_coeff(crv.coeff{1}, tind0/brkpnt(1));
		else
			pnt = evaluate_coeff(crv.coeff{2}, (tind0-brkpnt(1))/brkpnt(2));
		end
	else
		if numel(crv.coeff) ~= 1, error('Should be a polynom!'); end
		pnt = evaluate_coeff(crv.coeff{1}, tind);
	end

	return;
end
