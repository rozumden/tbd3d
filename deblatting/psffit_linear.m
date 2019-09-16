function [psf coeffs] = psffit_linear(psf, thresh)
% fis pw linear or quadratic polynomial to PSF and re-renders. Different initialization approach than psffit

% thresholding, find relevant pixel coordinates
bw = psfsupp(psf, thresh);
bw = imfill(bw, 'holes'); % to get rid of potential problems with euaclidian path winding around holes

% use euclidian tree to get 'clean skeleton'
m = bwdist(~bw); m = double(max(m(:))); % max 'thickness' of the psf
bw2 = imclose(bw, strel('disk', ceil(m*.75), 0)); % pre-cleaning edges of the trace
bw2 = bwmorph(bw2, 'thin', Inf); % required to get sinlge pixels paths, otherwise the euclid tree wiggles
[x1 x2] = find(bw2);
[P approx_len] = maxMinPath([x1 x2]); % shortest path in the previous skeleton
bw2 = false(size(bw2)); bw2((x2(P)-1)*size(bw2,1)+x1(P)) = true; % hopefully rather clean 'skeleton' of the original psf trace
path = [x1(P) x2(P)];

% cluster skeleton into 'k' groups to get some averaging and remove noise in the trajectory, hopefully along the original psf
k = max(ceil(nnz(bw2)/16),2); % kepp approx 16px in each cluster
[x1 x2] = find(bw2);
[~, C] = kmeans([x1 x2], k, 'Start', path(round(linspace(1,size(path,1), k)),:));

% prepare different kinds of fits
if(k > 2)
	% % establish ordering of clusters "along the line" - form a graph and find spanning tree, hopefully a good approximation
	% dist = sqrt(sum((reshape(C,[],1,2)-reshape(C,1,[],2)).^2,3)); % distances between all pairs of points
	% tri = delaunay(C(:,1), C(:,2));
	% idx1 = [tri(:,1); tri(:,1); tri(:,2)]; idx2 = [tri(:,2); tri(:,3); tri(:,3)]; % all pairs of vertices connected by edges
	% idx = (idx2-1).*size(dist,1)+idx1; % sub2ind in dist matrix
	% A = zeros(size(dist)); A(idx) = dist(idx); % adjecency matrix resulting from the delaunay triangulation
	% A = max(A, A.'); % symmetrize, since all edges were taken only one-way
	% tree = minspantree(graph(A));
	% len = distances(tree);
	% [v1,v2] = find(len==max(len(:)));
	% P = shortestpath(tree, v1(1), v2(1)); % longest path (set of indices into C) traversible in the tree - assume that this is the natural ordering of cluster points along the curve
	P = maxMinPath(C);
	C = C(P, :);

	% calculate angles between individual steps
	d = C(2:end,:)-C(1:end-1,:); % displacement vectors in each step along the path
	d = d./sqrt(sum(d.^2,2)); % normalized to unit length
	a = sum(d(2:end,:).*d(1:end-1,:),2); % angles in individual steps (cos(.), 1=straight step, 0=90 degrees, -1=step back)

	% identify breaking points
	% condition: is local max in angles, is at least some angle (non-straight)
	bends = find(ispeak(a) & a < .92);

	% different curve kinds
	% linear
	fit_params(1) = struct('control_pts', C([1 end],:), 'num_curves', 1, 'dbg', 'lin1');
	% quadratic
	% fit_params(1) = struct('control_pts', [C(1,:); C(floor((size(C,1)+1)/2),:); C(end,:)], 'num_curves', 1, 'dbg', 'quad1');

	% pw linear with 1 breaking point
	% (note: does not use 'bends' - is more stable when accidentally no bends are detected)
	% take first+last point as p0,p2 and calculate distance to the connecting line - take the point with gratest dist as p1
	p0 = C(1,:); p2 = C(end,:);
	v = p2-p0; v = v/sqrt(v*v');
	p = (C-p0)*[-v(2);v(1)]; % perpendicular projection
	[~,idx] = max(abs(p));
	p1 = C(idx,:);
	% fit_params(end+1) = struct('control_pts', [p0;p1;p2], 'num_curves', 2, 'dbg', 'lin2');

	if(numel(bends) >= 1)
		% pw quadratic with 1 breaking point
		[~,idx] = min(a(bends)); idx = bends(idx); % biggest bend (index into 'a')
		% quickfix: - do not allow bends in the outermost segments
		if(idx < 1 && idx < length(a))
			mid1 = C(floor((1+(idx+1))/2),:); % first midle control pt
			mid2 = C(floor((size(C,1)+(idx+1))/2),:); % second midle control pt
			% fit_params(end+1) = struct('control_pts', [C(1,:);mid1;C(idx+1,:);mid2;C(end,:)], 'num_curves', 2, 'dbg', 'quad2');
		end
		% if(idx > 1)
		% 	mid1 = C(floor((1+(idx+1))/2),:); % first midle control pt
		% else % bend in the first junction
		% 	mid1 = mean(C([1 2],:),1); % first midle control pt
		% end
		% if(idx < length(a))
		% 	mid2 = C(floor((size(C,1)+(idx+1))/2),:); % second midle control pt
		% else % bend in the last junction
		% 	mid2 = mean(C([end-1 end],:),1); % second midle control pt
		% end
		%fit_params(end+1) = struct('control_pts', [C(1,:);mid1;C(idx+1,:);mid2;C(end,:)], 'num_curves', 2, 'dbg', 'quad2');

		% pw linear with 2 breaking pts
		if(numel(bends) >= 2)
			[~,p] = sort(a(bends), 'ascend'); % take biggest 2 bends, bends(p(1:2))
			% fit_params(end+1) = struct('control_pts', C([1; sort(bends(p(1:2)))+1; end],:), 'num_curves', 3, 'dbg', 'lin3');
		end
	end
else % psf too short, try quadratic/linear fit and be done with it
	% linear
	fit_params(1) = struct('control_pts', C, 'num_curves', 1, 'dbg', 'lin1');
	% quadratic
	% fit_params(1) = struct('control_pts', [C(1,:); mean(C,1); C(2,:)], 'num_curves', 1, 'dbg', 'quad1');
end

% perform all fits
[x1 x2] = find(bw); x = [x1 x2];
for i=1:length(fit_params)
	sampling_rate = ceil(1*(approx_len/fit_params(i).num_curves));
	[fits(i).coeffs,~,err] = fitPwPolynomial(x, fit_params(i).control_pts, fit_params(i).num_curves, psf(bw), 20*fit_params(i).num_curves, 1e-2, sampling_rate);
	fits(i).err = mean(err(:));
	
	% % DBG - render
	% fits(i).img = trajRender(size(psf), cellfun(@(x)fliplr(x).', fits(i).coeffs, 'UniformOutput', false));
end

% select bes fit, render
[~,idx] = min([fits.err]);
psf = trajRender(size(psf), cellfun(@(x)fliplr(x).', fits(idx).coeffs, 'UniformOutput', false));
psf = psf/sum(psf(:));
coeffs = fits(idx).coeffs;


% % DBG
% disp(a.');
% % show all posibilities - control pts or fits
% for i=1:length(fit_params)
% 	figure(i);
% 	% imshow(psf,[]); title(fit_params(i).dbg);
% 	imshow(imoverlay(psf/max(psf(:)), fits(i).img, [1 0 0], .4)); title(sprintf('%s, err=%.2e', fit_params(i).dbg, fits(i).err));
% 	hold on;
% 	plot(C(:,2), C(:,1), '-b+');
% 	% control pts
% 	plot(fit_params(i).control_pts(:,2), fit_params(i).control_pts(:,1), '+g');
% 	% fit
% 	hold off;
% end

end

function [P len] = maxMinPath(x)
% finds min tree in a graph and longest traversable path in the tree - approximates traversing the curve
% x = [x1(:) x2(:)] coordinates

dist = sqrt(sum((reshape(x,[],1,2)-reshape(x,1,[],2)).^2,3)); % distances between all pairs of points
if(size(x,1) > 10) % not too few points
	try
		tri = delaunay(x(:,1), x(:,2));
		idx1 = [tri(:,1); tri(:,1); tri(:,2)]; idx2 = [tri(:,2); tri(:,3); tri(:,3)]; % all pairs of vertices connected by edges
		idx = (idx2-1).*size(dist,1)+idx1; % sub2ind in dist matrix
		A = zeros(size(dist)); A(idx) = dist(idx); % adjecency matrix resulting from the delaunay triangulation
		dist = max(A, A.'); % symmetrize, since all edges were taken only one-way
	catch
		% use dist matrix directly
	end
end
tree = minspantree(graph(dist));
len = distances(tree);
[v1,v2] = find(len==max(len(:)));
[P len] = shortestpath(tree, v1(1), v2(1));
end