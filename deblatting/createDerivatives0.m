function [Dx Dy] = createDerivatives0(sz2)
% returns sparse matrices that work as simple [-1 1] derivative operators
% does zero padding at the boundary
%
% sz2 - image 2-size

N = prod(sz2);
idx = reshape(1:N, sz2);
one = ones(N,1);

% x (forward differences with zero padding at boundary)
idx1 = idx(:, 1:end-1);
idx2 = idx(:, 2:end);
Dx = sparse([idx(:); idx1(:)], [idx(:); idx2(:)], [one; -one(1:numel(idx2))], N, N);

% y (forward differences with zero padding at boundary)
idx1 = idx(1:end-1, :);
idx2 = idx(2:end, :);
Dy = sparse([idx(:); idx1(:)], [idx(:); idx2(:)], [one; -one(1:numel(idx2))], N, N);
end
