function pnt = find_intersect(coeff0, coeff1, fit_iv, im)
Fit0 = myTrajRender(size2(im), coeff0, [fit_iv(1) fit_iv(2)+1]);
Fit1 = myTrajRender(size2(im), coeff1, [fit_iv(1) fit_iv(2)+1]);
[y0,x0] = find(Fit0);
[y1,x1] = find(Fit1);

[IDX,D] = knnsearch([x0 y0],[x1 y1]);
minid = find(D == min(D));

pnt = [x1(minid(end)) y1(minid(end))];
