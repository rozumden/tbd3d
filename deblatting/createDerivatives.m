function [Dx Dy] = createDerivatives(img_size_or_mask)
% returns sparse matrices that work as simple [-1 1] derivative operators (which is also surprisingly fast)
% does not differentiate at boundary (of image or mask) - returns zero there
%
% img_size_or_mask - either image 2-size, or logical mask of pixels where derivatives are to be computed. Operators are then square matrices with appropriate size - either prod(size) or nnz(mask)

if(islogical(img_size_or_mask))
	sz = size2(img_size_or_mask,1:2);
	mask = img_size_or_mask;
	N = nnz(mask);
else
	sz = img_size_or_mask(1:2);
	mask = true(sz);
	N = prod(sz);
end
idx = zeros(sz);
idx(mask) = 1:N;

% x (forward differences, fixed zeros at positions where one of the neighboring pixels lies outside of the mask)
idx1 = idx(:, 1:end-1);
idx2 = idx(:, 2:end);
keep = idx1 > 0 & idx2 > 0; idx1 = idx1(keep); idx2 = idx2(keep);
Dx = sparse([idx1; idx1], [idx2; idx1], [ones(numel(idx1),1); -ones(numel(idx1),1)], N, N);

% y (forward differences, fixed zeros at positions where one of the neighboring pixels lies outside of the mask)
idx1 = idx(1:end-1, :);
idx2 = idx(2:end, :);
keep = idx1 > 0 & idx2 > 0; idx1 = idx1(keep); idx2 = idx2(keep);
Dy = sparse([idx1; idx1], [idx2; idx1], [ones(numel(idx1),1); -ones(numel(idx1),1)], N, N);
end

% ORIGINAL VERSION (non-circular)
% function [Dx Dy] = createDerivatives(img_size_or_mask, matrix_sz)
% % returns sparse matrices that work as simple [-1 1] derivative operators (which is also surprisingly fast)
% % does not differentiate at boundary (of image or mask) - returns zero there
% %
% % img_size_or_mask - either image 2-size, or logical mask of pixels where derivatives are to be computed. Operators then have appropriate size - either prod(size) or nnz(mask)
% % matrix_sz - szie of output matrix either 'valid' (num_cols is nnz(mask), num_rows is such that only non-trivially-zero pixels are returned) or 'full' - num_cols=num_rows=prod(img_size)
% % FIXME: ^^ not implemented, square matrices with minimal input (=nnz(mask))

% if(islogical(img_size_or_mask))
% 	sz = size2(img_size_or_mask,1:2);
% 	mask = img_size_or_mask;
% 	N = nnz(mask);
% else
% 	sz = img_size_or_mask(1:2);
% 	mask = true(sz);
% 	N = prod(sz);
% end
% idx = zeros(sz);
% idx(mask) = 1:N;
% %idx = reshape(1:prod(sz), sz(1:2));

% % x (forward differences, fixed zeros for the rightmost column in the result)
% idx1 = idx(:, 1:end-1);
% idx2 = idx(:, 2:end);
% keep = idx1 > 0 & idx2 > 0; idx1 = idx1(keep); idx2 = idx2(keep);
% Dx = sparse([idx1; idx1], [idx2; idx1], [ones(numel(idx1),1); -ones(numel(idx1),1)], N, N);

% % y (forward differences, fixed zeros for the bottommost row in the result)
% idx1 = idx(1:end-1, :);
% idx2 = idx(2:end, :);
% keep = idx1 > 0 & idx2 > 0; idx1 = idx1(keep); idx2 = idx2(keep);
% Dy = sparse([idx1; idx1], [idx2; idx1], [ones(numel(idx1),1); -ones(numel(idx1),1)], N, N);
% end