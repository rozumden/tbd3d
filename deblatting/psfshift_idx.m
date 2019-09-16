function [idx_large] = psfshift_idx(small, sz_large)
% variant of psfshift intended for repeated use in a loop (subsequent calls are faster)
% determines index pairing between 'small' and 'large' image such that when large=psfshift(small, sz_large) then large(idx_large)==small(mask_small) (ie also optionally includes masking in the 'small' img, pixels not in mask are ignored and not psfshifted)
%
% small - either logical mask of relevant pixels in 'small' image, or size of the small image
if(islogical(small)) % mask
	temp = zeros(size(small));
	temp(small) = 1:nnz(small);
else % size
	temp = reshape(1:prod(small),small);
end
temp = psfshift(temp, sz_large);
idx_large = find(temp); temp_idx = temp(idx_large);
pos = zeros(size(temp_idx)); pos(temp_idx) = 1:numel(temp_idx); % permutation of indices 1:n in temp_idx
idx_large = idx_large(pos);
end