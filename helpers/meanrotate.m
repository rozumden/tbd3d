function [f0] = meanrotate(f)
ff = [];
stp = 10;
for ang = stp:stp:(360-stp)
	ff = cat(4, ff, imrotate(f, ang, 'bilinear', 'crop'));
end

f0 = mean(ff, 4);