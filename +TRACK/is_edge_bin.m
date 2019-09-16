function [res] = is_edge_bin(bin_c)
res = sum(bin_c(1,:)) + sum(bin_c(end,:)) + sum(bin_c(:,end)) + sum(bin_c(:,1));
