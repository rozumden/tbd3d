function res = size2(A, dims)
% similar to  `size` but allows for arbitrary dimension specification like 1:3, [3 1 2] etc
% returns 2d size by default, useful for images

if(nargin < 2)
	dims = 1:2; % returns 2-size by default
end

sz = size(A);
res = ones(size(dims)); % will automatically contain '1' for trivial (high) dimensions
valid = dims <= ndims(A); % non-trivial dimensions
res(valid) = sz(dims(valid));
end