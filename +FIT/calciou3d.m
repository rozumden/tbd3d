function iou = calciou3d(p1, p2, r)
dists = sqrt( sum((p1 - p2).^2,2) );
dists(dists > 2*r) = 2*r(dists > 2*r);

h = (r - dists./2); % height of the cap
A = (1/3)*pi.*(h.^2).*(3.*r - h); % volume of sphere cap

I = 2*A; % volume of intersection
U = 2*(4/3) * pi * r.^3 - I;
iou = I ./ U;
