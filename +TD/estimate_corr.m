function [pnts1, pnts2] = estimate_corr(flow, M, I1, I2, mag_th, grad_th)
grad1 = imgradient(rgb2gray(I1));
grad2 = imgradient(rgb2gray(I2));
grad = max(grad1,grad2);
inl = flow.Magnitude >= mag_th & grad >= grad_th;
inl(~M) = 0;
[y,x] = find(inl);
pnts1 = double([x y]);
pnts2 = double([x+flow.Vx(inl) y+flow.Vy(inl)]);
% Fmatrix = estimateFundamentalMatrix(pnts1, pnts2,'Method','RANSAC','NumTrials',2000,'DistanceThreshold',1e-4);
% [R, C, inlierIdx] = helperEstimateRelativePose(pnts1, pnts2, KP);
% [E, inlierIdx] = estimateEssentialMatrix(pnts1, pnts2, KP);
