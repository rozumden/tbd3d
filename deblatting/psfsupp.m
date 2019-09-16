function [keep psf] = psfsupp(psf, levels)
%PSFSUPP2 Rougly (based on thresholding) determines 'support' (significant pixels) of the given PSF.
%
% levels - decreasing values of thresholds (relative to the maximum value, as if 'psf' is fully stretched [0,1] image).
% Pixels with value grater than levels(1) are base of the support. At each higher (i-th) level, pixels value greater than levels(i) and connected to the previous ((i-1)-th) level are added to the support.
%
% keep - mask of PSF support
% psf - input PSF with values ouisde the determined support removed and normed to sum to 1

if(nargin < 2)
	levels = [.15 .05]; % [.2 .05] or [.1 .05] possible
end

% scale PSF to [0,1]
m = max(psf(:));
% if(m <= 0)
% 	keep = false(size(psf));
% 	psf = zeros(size(psf));
% 	return;
% end
% psf = psf / m;

keep = false(size(psf));
for i=1:length(levels)
	c = psf > levels(i)*m & ~keep; % candidates for inclusion
	if(i == 1) % first iteration
		add = c;
	else
		add = is_connected(add, c); % add only parts connected to previous layer of current support estimate
	end
	keep = keep | add;
	if(~any(add(:)))
		break;
	end
end

psf(~keep) = 0;
if(nnz(psf))
	psf = psf/sum(psf(:));
end
end

function c = is_connected(base, new)
% returns components of 'new' that are connectd to 'base'

labels = bwlabel(new);
boundary = imdilate(base, ones(3)) & (~base);

% pick components which have pixels close to the 'base'
keep = unique(labels(new & boundary));
c = false(size(base));
for i=1:length(keep)
	c = c | (labels == keep(i));
end
end
