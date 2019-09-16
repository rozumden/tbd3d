function [coor] = change_coor(mbb, x, y)
wdir = mbb(:,2) - mbb(:,1);
wdir = wdir / norm(wdir);
hdir = mbb(:,3) - mbb(:,2);
hdir = hdir / norm(hdir);

coor = mbb(:,1) + (x'-1).*hdir + (y'-1).*wdir;