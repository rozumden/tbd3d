function [peaks maxes mins] = ispeak(x, type)
% returns logical vector the size of x(=vector) with true in locations of local extrema. For 'flat' peaks returns true only in the first peak pixel.
%
% type: 'max' = local maximum, 'min' = local minimum, 'any' (default) = local max or min.

if(nargin < 2)
	type = ''; % any
end

switch(type)
	case 'max'
		leftmax = true(size(x)); % preallocation to avoid having to copy 'x' of unwnown size (orientation)
		leftmax(2:end) = x(1:end-1) < x(2:end);
		rightmax = true(size(x));
		rightmax(1:end-1) = x(2:end) <= x(1:end-1);

		peaks = leftmax & rightmax;
		maxes = peaks;
		mins = [];
	case 'min'
		leftmin = true(size(x));
		leftmin(2:end) = x(1:end-1) > x(2:end);
		rightmin = true(size(x));
		rightmin(1:end-1) = x(2:end) >= x(1:end-1);

		peaks = leftmin & rightmin;
		maxes = [];
		mins = peaks;
	otherwise
		leftmax = true(size(x));
		leftmax(2:end) = x(1:end-1) < x(2:end);
		rightmax = true(size(x));
		rightmax(1:end-1) = x(2:end) <= x(1:end-1);
		leftmin = true(size(x));
		leftmin(2:end) = x(1:end-1) > x(2:end);
		rightmin = true(size(x));
		rightmin(1:end-1) = x(2:end) >= x(1:end-1);

		maxes = leftmax & rightmax;
		mins = leftmin & rightmin;
		peaks = maxes | mins;
end
end
