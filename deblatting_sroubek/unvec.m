function image = unvec(v,y,x);
%
% image = unvec(v,y,x);
% Unvectorize an vector v of size (1 , y*x) and return matrix of size (x,y)
%

if nargin == 3
  s = [y x];
else
  s = y;
end
image = reshape(v,s);
