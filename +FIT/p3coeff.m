function [coeff] = p3coeff(p0,p,p1,t0)
if nargin < 4
	t0 = 0;
end

A1 = [1 0 t0 0; 0 1 0 t0; 1 0 t0+1 0; 0 1 0 t0+1]; 
b1 = [p0; p];
x1 = A1 \ b1;

A2 = [1 0 t0 0; 0 1 0 t0; 1 0 t0+1 0; 0 1 0 t0+1]; 
b2 = [p; p1];
x2 = A2 \ b2;

coeff = {[x1(1:2) x1(3:4)], [x2(1:2) x2(3:4)]};