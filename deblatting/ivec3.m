function [img] = ivec3(vec, sz)
%IVEC3 Utility, vectorized RGB image back to image (inverse of vec3)
img = reshape(vec, sz(1), sz(2), [], IF(numel(sz)>=4, @()sz(4),1)); % works with multiple stacked RGB images (return image with dims H-W-C-N)
end

