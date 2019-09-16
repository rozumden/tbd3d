function [srctable,xx,yy] = estimate_speeds(curves, r0, add_fctr)
if nargin < 3
	add_fctr = '';
end
srctable = '';
xx = []; yy = [];
for ci = 1:numel(curves)
	crv = curves(ci);
	if strcmp(crv.type,'joint') || strcmp(crv.type,'prediction')

		if numel(crv.coeff) > 1, error('Too many coeff'); end

		cf = crv.coeff{1};
		max_power = (numel(cf) / 2) - 1;

		addplotB = sprintf('\\addplot [domain=%d:%d,color=purple]',crv.fit_iv(1),crv.fit_iv(2)+1);

		eq1 = num2str(cf(1,2));
		eq2 = num2str(cf(2,2));

		for powi = 2:max_power
		    eq1 = [eq1 sprintf(' + %d*(%.12f)*x^%d',powi,cf(1,powi+1),powi-1)];
		    eq2 = [eq2 sprintf(' + %d*(%.12f)*x^%d',powi,cf(2,powi+1),powi-1)];
		end

		srctable = [srctable addplotB sprintf('({x%s},{(((%s)^2 + (%s)^2)^0.5)/%.4f}); \n',add_fctr,eq1,eq2,r0)];
		xx0 = crv.fit_iv(1):0.2:(crv.fit_iv(2)+1);
		xx = [xx xx0];
		for x = xx0
			yy = [yy norm(evaluate_dcoeff(cf, x))./r0];
		end
	elseif strcmp(crv.type,'bounce')
		crv0 = curves(ci-1);
		crv1 = curves(ci+1);
		p0 = norm(evaluate_dcoeff(crv0.coeff{1},crv0.fit_iv(2)+1));
		p1 = norm(evaluate_dcoeff(crv1.coeff{1},crv1.fit_iv(1)));
		t0 = crv.fit_iv(1);
		t1 = crv.fit_iv(2)+1;
		cf = [1 t0; 1 t1] \ [p0; p1];

		addplotB = sprintf('\\addplot [domain=%d:%d,color=purple]',t0,t1);
		srctable = [srctable addplotB sprintf('({x%s},{(%.4f + %.4f*x)/%.4f}); \n',add_fctr,cf(1), cf(2),r0)];
	
		xx = [xx linspace(t0,t1,10)];
		yy = [yy linspace(p0,p1,10)/r0];
	elseif strcmp(crv.type,'connect')
		disp('Skipping connect');
	else
		error('Unknown type');
	end
end
