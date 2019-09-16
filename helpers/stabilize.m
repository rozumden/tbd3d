function [im2p,tform] = stabilize(im1,im2)
im1_intensity = rgb2gray(im1);
im2_intensity = rgb2gray(im2);
ptThresh = 0.1;
pointsA = detectFASTFeatures(im1_intensity, 'MinContrast', ptThresh);
pointsB = detectFASTFeatures(im2_intensity, 'MinContrast', ptThresh);
	
% Extract FREAK descriptors for the corners
[featuresA, pointsA] = extractFeatures(im1_intensity, pointsA);
[featuresB, pointsB] = extractFeatures(im2_intensity, pointsB);

indexPairs = matchFeatures(featuresA, featuresB);
pointsA = pointsA(indexPairs(:, 1), :);
pointsB = pointsB(indexPairs(:, 2), :);

try 
	tform = estimateGeometricTransform(pointsB, pointsA, 'affine');
catch
	im2p = im2;
	return
end
im2p = imwarp(im2, tform, 'OutputView', imref2d(size(im2)));
