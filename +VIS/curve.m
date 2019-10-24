function [] = curve(curves)
if iscell(curves)
	step_size = 0.1;
	for ci = 1:(numel(curves))
		for cci = 1:numel(curves{ci})
			for tind = 0:step_size:1
				pnt = evaluate_coeff(curves{ci}{cci}, tind);
				plot(pnt(1),pnt(2),'.y');
			end
		end
	end	
else
	step_size = 0.05;
	tmin = curves(1).fit_iv(1);
	tmax = curves(end).fit_iv(2) + 1;
	last_pnt = [];
	for tind = tmin:step_size:tmax
		pnt = evaluate_curve_coeff(curves, tind);
		if isempty(last_pnt)
			plot(pnt(1),pnt(2),'.m');
		else
			plot([pnt(1) last_pnt(1)],[pnt(2) last_pnt(2)],'m','LineWidth',4);
		end
		last_pnt = pnt;
	end
end