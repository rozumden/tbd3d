function [pnt] = evaluate_dcoeff(coeff, ind)
max_power = size(coeff,2)-2;
t = ind.^[0:max_power];
cf = coeff(:,2:end).*[1:max_power+1];
pnt = t*cf.';
