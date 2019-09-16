function models = seqransac_lin(pts, w, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2, params)

% params
eps = 1e-6;
inlier_thresh = 1.5+eps;
% gap_thresh_small = 1.5+eps;
% gap_thresh_big = 5+eps;
% gap_thresh_nonsig = 2+eps; % gap threshold when looking for non-significant runs (2nd stage)
% sig_count1 = 4;
% sig_pcent1 = .4; % ^gamma
% sig_count2 = 6;
% sig_pcent2 = .8; % ^gamma
max_samples = params.seqransac_lin_max_samples; % 120pts full

% note: (HACK: removed and only sig segments are found, uncomment #alpha block to serach for non-sig segments in the second phase) first finds significant segments only (ie segment which starts and ends with is consecutive run of significant pts) and where no significant segments exist, finds remaining non-sig (with possibly stricter threshold for gap)

free_idx = 1:size(pts,1); % pts not yet in existing models (idx to pts)
% lbl_sum = accumarray(labels, w, [length(labels) 1]); % value of individual components
model_cs = {}; model_val = [];
sz = max(pts,[],1)+2; % psf size (cropped)
while(numel(free_idx) > 1)
	% look for currently best model
	num_samples = numel(free_idx)*(numel(free_idx)-1)/2; % total max number of samples available
	if(num_samples > max_samples)
		idx = randperm(num_samples, max_samples);
	else
		idx = 1:num_samples; % DBG - deterministic 'try all' approach
	end
	idx = nthtuple(idx, 2); % idx to free_idx or pts_free - pairs of pts for mss
	pts_free = pts(free_idx, :); w_free = w(free_idx); % labels_free = labels(free_idx);
	new_val = 0;

	% recalculate component labeling (connectivity between remaining pieces of psf)
	temp = accumarray(pts_free, 1, sz); temp = temp | imclose(temp, ones(3)); % psf super-mask of the remaining pts
	[labels num] = bwlabel(temp); labels = labels((pts_free(:,2)-1)*sz(1)+pts_free(:,1)); % component label of each remaining pt
	lbl_sum = accumarray(labels, w_free, [num 1]); % value of individual components

	for i = 1:size(idx,1)
		% distance to line - general consensus set
		mss = pts_free(idx(i,:), :);
		d = (mss(2,:)-mss(1,:)); d = d/sqrt(d*d(:)); % unit line direction vector
		cs = find(abs((pts_free-mss(1,:))*[-d(2); d(1)]) <= inlier_thresh); % concensus set (potential model inliers); distance of pts from the line below threshold
		t = (pts_free(cs,:)-mss(1,:))*d(:); % distance parametrization of all cs pts on the mss line

		% find best run
		[cs2 val] = findSignificantRun(t, w_free(cs), labels(cs), lbl_sum, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2);
		if(val > new_val)
			new_val = val; new_cs = cs(cs2); % currently best model
		end
	end
	if(new_val == 0) % no non-single-px model could be found
% 		if(sig_thresh) % all significant runs exhausted, find remaining non-sig runs #alpha
% 			sig_thresh = 0;
% 			gap_thresh_small = gap_thresh_nonsig;
% 			continue;
% 		end
		break;
	end 
	temp = free_idx(new_cs); free_idx(new_cs) = []; new_cs = temp; % recalculate to global idx to pts and remove from free

	% add new model to existing collection
	model_cs{end+1} = new_cs; model_val(end+1) = new_val;
	
% 	% DBG
% 	figure(1);
% 	imshow(psf); hold on;
% 	for i=1:length(model_cs)
% 		plot(pts(model_cs{i},2), pts(model_cs{i},1), 'o');
% 	end
% 	hold off;
% 	pause
end

% refit models, ~~return sorted best->worst
% [~, p] = sort(model_val, 'descend'); note: sorting removed
models = struct([]);
for i=1:length(model_cs)
	idx = i; % idx = p(i); note: sorting removed
	[end_pts d] = fitsegment2(pts(model_cs{idx},:)); % refit to get segment params
	
	models(i).cs = model_cs{idx}; % consensus set (idx to pts)
	models(i).val = model_val(idx); % value
	models(i).d = d; % unit direction vector
	models(i).end_pts = end_pts; % end points (start;end)
	models(i).len = sqrt(sum((end_pts(2,:)-end_pts(1,:)).^2)); % segment length
end
end

% function [val cs] = find_cs(mss, pts, w, inlier_thresh, gap_thresh)
% % find consensus set and correspoinding 'model' (particular line segment) to the given minimal sample set (2 pts defining a line)

% % distance to line - general consensus set
% d = (mss(2,:)-mss(1,:)); d = d/sqrt(d*d(:)); % unit line direction vector
% %factor = max(abs(d)); % normalizing factor to treat distances between pts in a diagonal and axis-aligned line equally
% factor = 1; % seems to work better than ^^
% cs = find(abs((pts-mss(1,:))*[-d(2); d(1)]) <= inlier_thresh/factor); % concensus set (potential model inliers); distance of pts from the line below threshold

% % project all pts onto the line, find best uninterrrupted run
% t = (pts(cs,:)-mss(1,:))*d(:); % distance parametrization of all cs pts on the mss line

% [cs val] = findSignificantRun(t, w(cs), is_sig(cs), gap_thresh_small, gap_thresh_big, min_sig_count)


% [t, p] = sort(t); cs = cs(p);
% dt = t(2:end)-t(1:end-1); % spacing between i+1-th and i-th pt in cs
% gaps = [0; find(dt > gap_thresh/factor); numel(t)]; % i-th consecutive run is gaps(i)+1:gaps(i+1)
% cum_val = [0; cumsum(w(cs))]; % cumulative value of pts along the line
% vals = cum_val(gaps(2:end)+1)-cum_val(gaps(1:end-1)+1); % values of consecutive runs
% num = gaps(2:end)-gaps(1:end-1); % number of pts in a run
% [val ix] = max(vals.*(num > 1)); % run with max value (ignore single-pixel "runs")
% cs = cs(gaps(ix)+1:gaps(ix+1)); % idx to pt, consensus set
% end

function [end_pts d t] = fitsegment2(pts)
% note: very similar results to fitsegment, significant difference only in near-symmetrical cases (points on square etc) (this version is slightly faster than ~1)
% fit line by pca
m = mean(pts,1); pts0 = pts - m;
c = (pts0.'*pts0); % covariance
[V,D] = eig(c);
[~, p] = sort(diag(D));
d = V(:, p(2)).'; % largest eigenvector, unit direction vector
t = pts0*d(:); [t1 t2] = bounds(t); % projections onto the principal direction - min and max (furthest two pts)
end_pts = [m+t1*d; m+t2*d];
end
