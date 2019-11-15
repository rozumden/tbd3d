function [] = curve3d(curves, szs, ind)
step_size = 0.5;
tmin = curves(1).fit_iv(1);
tmax = curves(end).fit_iv(2) + 1;
last_pnt = [];
for tind = tmin:step_size:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	pnt = [pnt interp1(ind,szs,tind)];
	if isempty(last_pnt)
		plot3(pnt(1),pnt(2),pnt(3),'.g'); hold on
	else
		plot3([pnt(1) last_pnt(1)],[pnt(2) last_pnt(2)],[pnt(3) last_pnt(3)],'g','LineWidth',4);
	end
	last_pnt = pnt;
end
