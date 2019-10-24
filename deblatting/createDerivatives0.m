function [Dx Dy] = createDerivatives0(sz)
% Returns sparse matrices that work as simple [1 -1] derivative operators with zero padding at all boundaries.
% The actual derivative is one-pixel larger than the input image (essentially calculates deriavtives between pixels, incl. first and last) and is zero-padded (at bottom and right) so that Dx and Dy have same height.
% Both Dx and Dy are rectangular (tall) matrices.
%
% sz - image 2-size

N_in = prod(sz(1:2)); % input size, matrix width
N_out = prod(sz(1:2)+1); % output size, matrix height due to input and output zero padding
idx_in = reshape(1:N_in, sz(1:2));
idx_out = reshape(1:N_out, sz(1:2)+1);

% x - difference in the middle, one-sided 'difference' (pixel value against zero) at horizontal boundary; output zero padding at the bottom
v1 = ones(sz(1), sz(2)-1); v2 = ones(sz(1),1);
Dx = sparse(idx_out(1:end-1,[1 2:end-1 2:end-1 end]), idx_in(:,[1 1:end-1 2:end end]), [v2 -v1 v1 -v2], N_out, N_in);

% y - difference in the middle, one-sided 'difference' (pixel value against zero) at vertical boundary; output zero padding on the right
v1 = ones(sz(1)-1, sz(2)); v2 = ones(1,sz(2));
Dy = sparse(idx_out([1 2:end-1 2:end-1 end],1:end-1), idx_in([1 1:end-1 2:end end],:), [v2; -v1; v1; -v2], N_out, N_in);
end
