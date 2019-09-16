function h = psfshift(h, usize)
%PSFSHIFT Moves PSF center to origin and extends the PSF to be the same size as image (for use with FT). ipsfshift does the reverse.
%
% h				PSF to be shifted
% usize			target size (after extending with zeros)

% In conv2(., h, 'same'), the center pixel of 'h' has coordinates ceil((size(h)+1)/2), eg [3 3] for 4x4 array. In FT the center pixel has coords [1 1] (in MATLAB).

hsize = size(h); hsize = hsize(1:2);
h = circshift(padarray(h, usize(1:2)-hsize, 0, 'post'), -ceil((hsize+1)/2)+1);
end