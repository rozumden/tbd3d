function img = renderpwpoly2(coeffs, sz)
% fast rendering of piecewise linear or quadratic polynomial curves
% (could be extended to cubics, but higher order requires solving > quadratic eq for determining the splitting pts)
% 
% coeffs - cell array of coeffs for each curve. Each coeffs{i} is Nx2 array where N=2 for linear curves, N=3 for quadratic, inceasing order top to bottom

img = zeros(sz);

% algorithm - split each qudaratic curve into sufficiently many linear segments, which can be rendered easily
for i=1:length(coeffs)
	splits = render_poly2(coeffs{i});
	
	% render line segments between splits
	pts = round(splits.^(0:size(coeffs{i},1)-1)*coeffs{i}); % pts on the rendered polyline
	count = max(max(abs(pts(2:end,:)-pts(1:end-1,:)), [], 2)+1,2); % number of discrete rendered pts on each line segment (at least 2 so that it renders single pts as well)
	t = (0:max(count)-1).'; % preallocated timesteps
	for j=1:length(count)
		idx = round(pts(j,:) + (t(1:count(j))/(count(j)-1)).*(pts(j+1,:)-pts(j,:))); % sample line segment
		idx = idx(all(idx > 0 & idx <= sz,2), :); % discard indices outside the canvas size
		img((idx(:,2)-1)*sz(1) + idx(:,1)) = 1;
	end
end
end

function splits = render_poly2(coeffs)
splits = [0; 1]; % parametrization of split pts

% split quadratic curve
N = size(coeffs,1)-1; % poly order
if(N > 1)
	head_idx = 1; % index of the examined segment. head_idx == i indicates the curve segment between splits(i) and splits(i+1) must be examind while all splits before can be approxiamted linearly

	% evaluation at first split pt
	a1 = splits(head_idx).^(0:N)*coeffs;

	while(head_idx < length(splits))
		% evaluate curve at current split pts
		a2 = splits(head_idx+1).^(0:size(coeffs,2))*coeffs; % current next split
		d=a2-a1; % direction vector between split pts
		
		% determin single potential split pt as the furthes point from the approximatng line. Note: this is tailoer for quadratic curves
		val1 = coeffs(2,2)*d(1)-coeffs(2,1)*d(2); % appears in the eq for furthets point and determines of the line is geometrically quadratic (numerator)
		val2 = 2*(coeffs(3,1)*d(2)-coeffs(3,2)*d(1)); % appears in the eq for furthets point and determines of the line is geometrically quadratic (denominator)
		if(abs(val2) > abs(val1))
			t = val1/val2; % parametrization of the furthes pt
			% determine distance from the straight line
			p = t.^(0:N)*coeffs; % furthest pt, potential new split
			if(abs((p-a1)*[d(2); -d(1)]/sqrt(sum(d.^2))) > .5) % distance from line is above tolerance, split
				splits = [splits(1:head_idx); t; splits(head_idx+1:end)]; % add new split pt, keep head as is
				continue;
			end
		end

		% advance
		head_idx = head_idx + 1;
		a1 = splits(head_idx).^(0:N)*coeffs; % reevaluate
	end
end
end