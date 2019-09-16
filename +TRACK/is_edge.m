function [res] = is_edge(im_c, bgr_c)
delta = abs(im_c - bgr_c);
bin_c = sum(delta,3) > 9/255;
sumall = sum(bin_c(1,:)) + sum(bin_c(end,:)) + sum(bin_c(:,end)) + sum(bin_c(:,1));
res = sumall > 0.1*min(size(bin_c));