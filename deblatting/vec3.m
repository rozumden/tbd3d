function vec = vec3(img)
%VEC3 Utility, vectorize RGB image into columns, one per each color channel (like (:) but trets channels separately). The inverse is ivec3.
vec = reshape(img, [], size(img,3)*size(img,4)); % works for multiple stacked RGB images as well (stacking then works as [RGB1 RGB2 RGB3...])
end
