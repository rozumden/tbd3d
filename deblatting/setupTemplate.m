function [img m0 lambda_m0] = setupTemplate(img, pad_pcent, lambda_m0_in, lambda_m0_out)
if(nargin < 2) pad_pcent = 0; end
if(nargin < 3) lambda_m0_in = 0; end
if(nargin < 4) lambda_m0_out = 0; end

sz0 = size2(img);
pad = floor(sz0*pad_pcent);
sz = sz0 + 2*pad;

m0 = diskMask(sz, min(sz0)/2); % inner disk mask
img = padarray(imgaussfilt(img, max(min(size2(img))/18,1)), pad, 0, 'both'); % +blur
%lambda_m0 = zeros(sz); lambda_m0(m0) = lambda_m0_in; lambda_m0(~m0) = lambda_m0_out;
lambda_m0 = padarray(lambda_m0_in*ones(sz0), pad, lambda_m0_out, 'both');
end