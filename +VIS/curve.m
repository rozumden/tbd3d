function [] = curve(curves)
tmin = curves(1).fit_iv(1);
tmax = curves(end).fit_iv(2) + 1;
for tind = tmin:0.1:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	plot(pnt(1),pnt(2),'.r');
end
