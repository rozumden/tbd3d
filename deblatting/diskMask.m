function mask = diskMask(sz, r)
% disk mask of of size 'sz' (2x1) and radius 'r' (scalar). If one is not specified, then r=min(sz)/2 and vice versa.
% mask with 2*r = sz is exact fit
% perhaps counterintuitively returns even-sized masks for integer 'r', but the result is more realistic

% determine size/r (note: r is exactly size/2 (up to rounding), therefore for larger sizes, there is no singular 'dot' on the horizontal and vertical axes of the disk, which appears with the typical size=(2*r+1))
if(nargin < 2 || isempty(r))
	sz = sz([1 min(2,end)]);
	r = min(sz(1:2))/2;
elseif(isempty(sz))
	sz = round(2*r([1 1])); % size wrt 'r' rounded to .5
end

[x y] = meshgrid(-(sz(2)-1)/2:(sz(2)-1)/2, -(sz(1)-1)/2:(sz(1)-1)/2); % center pixel is (0,0) or no center pixel and 0.5 steps. 'r' is from center to the side (not center of the side pixels)
mask = x.^2 + y.^2 <= r^2;
end