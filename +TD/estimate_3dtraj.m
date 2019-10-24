function [szs,matF,matM,ind] = estimate_3dtraj(video, curves, M, n, params_tbd3d)
if ~exist('params_tbd3d','var')
	params_tbd3d.do_hier = true;
	params_tbd3d.iter_smoothing = 4;
end
params_tbd3d.do_intervals = IF(isfield(params_tbd3d, 'do_intervals'), @()params_tbd3d.do_intervals, false); 

F = ones([size(M) 3]);
[f,m] = estimateFM_mc(video, curves, M, F); 
% m = m.*M; f = f.*M;

if iscell(curves) %% GT
	if params_tbd3d.do_hier
		[matF, matM, ind] = TD.get_views_gt(video, curves, m, f, n, params_tbd3d);
	else
		[matF, matM, ind] = TD.get_views_gt_nonhier(video, curves, m, f, n, params_tbd3d);
	end
else
	[matF, matM, ind] = TD.get_views_curves(video, curves, m, f, params_tbd3d);
end
if size(matM,4) == 1
	matM = permute(repmat(matM,[1 1 1 1]),[1 2 4 3]);
end
% matF = matF.*M; matM = matM.*M;
[sz] = TD.estimate_3d(matF,matM,2);

if params_tbd3d.iter_smoothing > 0
	szs = smooth(sz, 'rlowess');
	for loopk = 1:params_tbd3d.iter_smoothing
		szs = smooth(szs, 'rlowess');
	end
end

% plot(ind, szs);

if ~params_tbd3d.do_intervals
	szs = reshape(szs, n, []);
	szs = num2cell(szs,1);
end

[ind, unisort] = unique(ind);
szs = szs(unisort);
matF = matF(:,:,:,unisort);
matM = matM(:,:,:,unisort);
