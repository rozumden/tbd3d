function iou = calciou(p1, p2, r)
dists = sqrt( sum((p1 - p2).^2,2) );
dists(dists > 2*r) = 2*r;

theta = 2*acos( dists./ (2*r) );
A = (r^2/2) * (theta - sin(theta)); % area of circular segment
I = 2*A; % area of intersection
U = 2* pi * r^2 - I;
iou = I ./ U;
