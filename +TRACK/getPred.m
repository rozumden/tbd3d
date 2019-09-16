function [bb, pts] = getPred(coeff0, ext, expos, sz, direc)
pts = FIT.trajPredictor(coeff0, direc, expos);
p0 = round(min(pts) - ext);
p0(p0 < 1) = 1;
p1 = round(max(pts) + ext);
szf = fliplr(sz);
p1(p1 > szf) = szf(p1 > szf);
bb = [p0 p1];