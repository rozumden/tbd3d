function [] = curve(curves)
step_size = 0.1;
if iscell(curves)
	for ci = 1:(numel(curves)+1)
		for tind = 0:step_size:1
			pnt = evaluate_vcoeff(curves{ci}, tind);
			plot(pnt(1),pnt(2),'.y');
		end
	end	
else
	tmin = curves(1).fit_iv(1);
	tmax = curves(end).fit_iv(2) + 1;
	for tind = tmin:step_size:tmax
		pnt = evaluate_curve_coeff(curves, tind);
		plot(pnt(1),pnt(2),'.r');
	end
end