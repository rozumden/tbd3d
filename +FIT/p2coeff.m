function [coeff] = p2coeff(p0,p1,t0,t1)

A1 = [1 0 t0 0; 0 1 0 t0; 1 0 t1 0; 0 1 0 t1]; 
b1 = [p0; p1];
x1 = A1 \ b1;

coeff = {[x1(1:2) x1(3:4)]};