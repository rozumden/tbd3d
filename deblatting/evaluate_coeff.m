function [pnt] = evaluate_coeff(coeff, ind)
max_power = size(coeff,2)-1;
t = ind.^[0:max_power];
pnt = t*coeff.';
