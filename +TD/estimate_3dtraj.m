function [szs,matF,matM] = estimate_3dtraj(video, curves, M)
F = ones([size(M) 3]);
[f,m] = estimateFM_mc(video, curves, M, F);
m = m.*M; f = f.*M;
[matF, matM, ind] = TD.get_views_curves(video, curves, m, f, true);
matF = matF.*M; matM = matM.*M;
[sz] = TD.estimate_3d(matF,matM,2);

szs = smooth(sz, 'rlowess');
plot(ind, szs);