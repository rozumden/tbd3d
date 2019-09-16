function [T, coeff, s, len] = getpoint(h, hext, th, mval)
[y x] = find(hext == mval);
y = y(1); x = x(1);
s = mval;
len = 1;
coeff = {[y x; 0.5 0.5]};
T = myTrajRender(size(h), cellfun(@(x)fliplr(x).', coeff, 'UniformOutput', false));
T = T / sum(T(:));