function [im_c, bgr_c] = make_crop_axis_parallel(im, bgr, bb)
im_c = im(bb(2):bb(4),bb(1):bb(3),:);
bgr_c = bgr(bb(2):bb(4),bb(1):bb(3),:);