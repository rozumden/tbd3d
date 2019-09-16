function [coor] = change_coor_back(mbb, xc, yc)
wdir = mbb(:,2) - mbb(:,1);
wdir = wdir / norm(wdir);
hdir = mbb(:,3) - mbb(:,2);
hdir = hdir / norm(hdir);

coor = [[xc(:)]'; [yc(:)]'] - mbb(:,1);
coor = [ coor'*wdir, coor'*hdir]+1;
