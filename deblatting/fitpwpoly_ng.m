function [coeffs info] = fitpwpoly_ng(pts, control_pts, poly_order, w, maxiter, sampling_delta, inlier_thresh, curvature_w)
% FIXME: udpatovat
% poly-order - scalar, 1-linear, 2-quadratic etc
%
% coeffs - poly coeffitients for each curve (cell array, coeffs{i} is matrix of coeffs for ith curve)
% info:
% dist2 - approx squared distance of input pt to curve
% curve_t - matching of input adn curve pts - curve_t(i) is parametrization of i-th input pt in its curve (see curve_idx)
% curve_idx - matching of input pts to curve - 1-based idx of closest curve for each input pt

% TODO:
% - vyradit body s w=0 na zacatku - spis napsat at se takove body nedavaji na vstup, at pak nemusim resit spravnou indexaci na vystupu (u dist2 apod)
% - zkusit nejak poznat a vyresit pripady kdy to z principu nefunguje (hlidat situace kdy krivka zdegeneruje do bodu apod)

% TODO: akutne

% aktualne (2/2019) nejlepsi (hlavne nejspolehlivejsi) verze ze vsech fitpwpoly* funkci
% NOTE: funguje dobre, krom zatim neimplementovanych veci (vracet vic info nez jen coeffs, hlidani singularity v prubehu)
% NOTE: znamy problem: nekdy se stane ze to 'dokonverguje' a zacne nepatrne oscilovat. Pravdepodobne to zpusobuje L2 corner loss (idx_boundary atd) tim ze bod preskakuje dobnitr a ven, ale uplne jsem to nezkoumal a neznam presnou pricinu. Zda se ze pomohlo kdyz se corner loss trochu posunul a krivka muze "prestelit" (viz #bravo) ale je mozne ze to je jen nahoda, navic i tak nekdy oscilace nastane, ale casteji to dokonvergovalo. Kazdopadne kdyz se dost nizko omezi pocet iteraci aby to nejelo vecnost tak by to v praxi snad nemelo vadit, protoze oscilace jsou nepatrne.
% NOTE: clen s omezovanim curvature (curvature_w apod) byl pridan dost narychlo a neni uplne otestovane. Kdp to neni striktne vzato 'curvature' (byl jsem liny to tak udelat) ale penalizuje to druhou derivaci (konkretne druha^2/prvni^4 - vychazi z orezaneho vztahu pro curvature ~ |druha|/|prvni|^2, chybi tam ta normala a je to invariantni jen k linearni reparametrizaci), coz neni uplne to same. Bohuzel je trochu zavisle na samplovani krivky (samplding_delta apod) protoze normovani neni uplne asi spravne takze vhodne curvature_w se trochu lisi podle tohoto nastaveni, ale neni to zas tolik a i tak je to v praxi pouzitelne

num_curves = (size(control_pts,1)-1)/poly_order;
abs_tol = 1e-3; % max between iteration control pt displacement for exit (originally 1e-3)
reset_tol = 2; % max between iteration control pt displacement to continue with current sampling setup (reset if diff>tol)
epsilon = .5; % irls
boundary_pad = .1; % distance; padding between actual curve end and first/last sampled point (same distance for all curves). This is important for correctly assinging input points to curves at corner and it is necessary that this padding is of const length across all curves (ie independed of curve length)
if(poly_order == 1)
	curvature_w = 0;
else
	curvature_w_coeffs = [1; -4; 4]; % polynomial coeffs that specify weighting based on curve parametrization (reaches 1 at t=0,1 and 0 at t=.5)
end

V = inv((linspace(0,1,poly_order+1).').^(0:poly_order)); % inverse of the vandermonde matrix that evaluates the polynomial at control pts, V transforms control pts to poly coefficients, therefore evaluation of polynomial at sampling pts is x = T*V*control_pts where T is the sampling vand matrix
control_idx = 1:poly_order:size(control_pts,1); % conrrespondence between curve_idx and its control pts: control_pts(control_idx(k):control_idx(k+1), :) are control pts of k-th curve

% preallocations
lhs = zeros(size(control_pts,1),size(control_pts,1)); % LHS matrix of the linear problem in each iter; contains subproblems for all curve segments with the eqs for shared points summed
rhs = zeros(size(control_pts)); % RHS for the main linear system. Again, eqs for all control points incl. sharing
idx_boundary = zeros(num_curves,2); % (i,:) is idx of (first,last) pt assigned to i-th curve (wrt ordering based on t-parametrization) - this is used for 'corner loss'(binding corner points to curve points); set to 0 when points exist on both sides of curve ends (the corner loss doesn apply) and to -1 is no curve points exist (this is used as a flag)
idx_curve = cell(1,num_curves); % indices of points beloging to i-th curve (inverse of curve_idx)
diff2 = Inf;

% main loop
for iter = 1:maxiter
	control_old = control_pts;

	% SETUP SAMPLING (reset sampling if a significant change occured)
	if(diff2 > reset_tol^2)
		% determine approximate curve length (note: this could reuse current sampled points but then it would always be 1 iteration old, which can be significant in early iterations)
		if(poly_order > 1)
			curve_len = zeros(num_curves,1); % curve_len(k) is approximate length of k-th curve
			% trace curves
			temp = ((linspace(0,1,poly_order*5).').^(0:poly_order))*V;
			for i=1:num_curves
				temp_pts = temp*control_pts(control_idx(i):control_idx(i+1),:);
				curve_len(i) = sum(sqrt(sum((temp_pts(2:end,:)-temp_pts(1:end-1,:)).^2,2)));
			end
		else % simple for pw lin curves
			curve_len = sqrt(sum((control_pts(2:end,:)-control_pts(1:end-1,:)).^2,2)); % curve_len(k) is length of k-th curve
		end

		% preallocations and indexing setup
		num_samples = max(ceil(curve_len./sampling_delta)+1,2); % # of sampled pts at each curve
		curve_pts = zeros(sum(num_samples), size(pts,2)); % preallocation - sampled pts of individual curves
		curve_dt = zeros(size(curve_pts,1), size(pts,2)); % preallocation - derivatives at sampled pts
		curve_pts_idx = [1; cumsum(num_samples)+1]; % indexing into curve_pts(.,:) for individual curves; curve_pts(curve_pts_idx(k):curve_pts_idx(k+1)-1, :) are sampled pts belonging to k-th curve
		is_first_last = false(size(curve_pts,1),3); is_first_last(curve_pts_idx(1:end-1),[1 3]) = 1; is_first_last(curve_pts_idx(2:end)-1,[2 3]) = 1; % indicator if given sampled pt is first or last in the curve. (i,:) is [is-fisrt is-last is-first-or-last] in i-th curve
		pt2curve_idx = cumsum(is_first_last(:,1)); % for curve point curve_pts(idx,:) the corresponding curve-idx is pt2curve_idx(idx)
		sampling_t = zeros(sum(num_samples),1);

		% determine offset of the sampling - slight padding between actual curve end (t=0,1) and first/last sampled point. This is important for correctly assinging input points to curves at corner and it is necessary that this padding is of const length across all curves (ie independed of curve length)
		sampling_offsets = zeros(num_curves,2); % (i,:) is (from,to) parameter values for i-th curve, corresponding to first and last sampled pt. Sth like (epsilon,1-epsilon)
		for i=1:num_curves
			sampling_offsets(i,:) = boundary_pad./sqrt(sum(([V(2,:); (1:poly_order)*V(2:end,:)]*control_pts(control_idx(i):control_idx(i+1),:)).^2,2)); % derivative sampled at curve ends, offset then specifies sampling 't' a given distance (boundary_pad) from the curve end
		end
		t_offset_min = -sampling_offsets.*(num_samples-1)./(1-sum(sampling_offsets,2)); % [start end] min values of t_offset allowed for each curve (to compensate for sampling_offsets and allow curve_t reaching but not exceeding 0 resp 1) (note: this formula assumes equi-parametric sampling, ie sampling_t=linspace(...) and needs to be updated if this is chnaged #alfa)

		for i=1:num_curves
			t = linspace(sampling_offsets(i,1), 1-sampling_offsets(i,2), num_samples(i)).'; % note: is this equi-parametric sampling is changed, the line marked #alfa needs to be updated
			sampling_t(curve_pts_idx(i):curve_pts_idx(i+1)-1) = t;
			sampling_V{i} = (t.^(0:poly_order))*V;
			sampling_V_dt{i} = ((t.^(0:poly_order-1)).*(1:poly_order))*V(2:end,:);
			if(curvature_w(1) > 0)
				dt = sampling_V_dt{i}*control_pts(control_idx(i):control_idx(i+1),:); % derivatives (FIXME: za chvili se bude samplovat znova, slo by si zapamatovat)
				d2t = ((t.^(0:poly_order-2).*(1:poly_order-1).*(2:poly_order))*V(3:end,:))./sum(dt.^2,2); % 2nd derivative at sampled pts, normed by dt^2 (the minimized quantity is second^2/first^4 per sample - this is before being squared) Note: this is not exactly geometric curvature, rather some simplification using second derivative only.
				if(numel(curvature_w) > 1) % optional increased weight at curve ends
					d2t = d2t.*(1+curvature_w(2)*(t.^(0:2)*curvature_w_coeffs).^2); % optional increased weighting at curve ends
				end
				curvature_term{i} = d2t.'*d2t/num_samples(i)^2;
			else
				curvature_term{i} = 0;
			end
		end
	end
	
	% trace curves
	for i=1:num_curves
		curve_pts(curve_pts_idx(i):curve_pts_idx(i+1)-1, :) = sampling_V{i}*control_pts(control_idx(i):control_idx(i+1),:);
		curve_dt(curve_pts_idx(i):curve_pts_idx(i+1)-1, :) = sampling_V_dt{i}*control_pts(control_idx(i):control_idx(i+1),:); % derivatives
	end

	% match each pt to a curve via shortest distance
	dist2 = sum((reshape(pts,size(pts,1),1,size(pts,2))-reshape(curve_pts,1,size(curve_pts,1), size(curve_pts,2))).^2,3); % (squared) distance "input_point-curve_point", all mutual combinations; indexed as (input_pts,curve_pts)
	[~,idx] = min(dist2,[],2); % closest curve point to each input point; idx needs to be recalculated to curve_idx/time_idx
	
	% assign pts to curves and vice versa (determine which pts belong to given curve)
	curve_idx = pt2curve_idx(idx); % idx of curve to which each input pt belongs
	% idx_curve = sparse(1:size(pts,1), curve_idx, true, size(pts,1),num_curves); % inverse of curve_idx; idx_curve(:,i) is logical array of pts belonging to i-th curve (this 'transposed' shape should fit better with matlab's storage of sparce matrices)
	for i=1:num_curves idx_curve{i} = find(curve_idx == i); end % faster than arrayfun etc and more convenient to use than logical mask via sparse(...)

	% find newighboring pt to idx-th pt so that I have two clsoest pt to each input pt and I can do projection on the linear appxorimation between the two pts
	p1 = curve_pts(idx,:);
	d = sum((pts-p1).*curve_dt(idx,:),2); % dot product with curve dt - determines direction in which the actually closest pt lies from the ampled closest pt
	idx2 = ones(size(idx)); idx2((d < 0 | is_first_last(idx,2)) & ~is_first_last(idx,1)) = -1; % offset to idx, index of the neighbor on the side where the projection lies: -1=left of idx, +1=right of idx
	
	% find the actual closest curve pt by linear approximation of the curve between 1st and 2nd determined closest pts
	p2 = curve_pts(idx+idx2, :); % line segment on which to project
	t1 = sampling_t(idx); t2 = sampling_t(idx+idx2); % timestaps of the nearest segment
	t_offset = sum((pts-p1).*(p2-p1),2)./sum((p2-p1).^2,2); % relative parametrization offset; t1+t_offset.*(t2-t1) should be the new parametrization (curve_t) of each pt (corresponds to orthogonal projection); requires clipping to [0,1]
	
	% keep original unclipped parametrization (used to reparametrize outermost curves to quickly match curve-ends exactly)
	t_offset_unclipped = t_offset;
	
	% clip offset to get actual curve_t
	t_offset(t_offset > 1) = 1; % upper bound - clip to '1' everywhere
	t_offset(~is_first_last(idx,3) & t_offset < 0) = 0; % lower bound - clip to 0 everywhere except at curve ends (both outermost ends and corners)
	temp = is_first_last(idx,1); t_offset(temp) = max(t_offset(temp), t_offset_min(curve_idx(temp),1)); % lower bound at beggining of curve - clip such that curve_t can reach t=0 (hence slightly negative values allowed for t_offset because of sampling_offset (t1 for 'first' curve point has t>0))
	temp = is_first_last(idx,2); t_offset(temp) = max(t_offset(temp), t_offset_min(curve_idx(temp),2)); % lower bound at end of curve - clip such that curve_t can reach t=1
	curve_t = t1 + t_offset.*(t2-t1); % best guess of the point->curve projection

	% recalculate dist2 quickly from t_offset - used for l1 reviewting
	dist2 = sum((pts - p2.*t_offset - p1.*(1-t_offset)).^2,2); % covnex comb of the two closest curve pts and distance there
	inliers = dist2 <= inlier_thresh^2; % (FIXME: the rest of the code will throw error if inliers is empty, maybe do some check...)
	
	% convergence - exit after pt matching so that additional info can be calculated
	if(diff2 < abs_tol^2)
		break;
	end

	% detect inliers and recalc all indexings to effectively discard outliers 
	for i=1:num_curves idx_curve{i} = idx_curve{i}(inliers(idx_curve{i})); end
	
	% find outermost curve ends for reparametrization
	temp = min(idx(inliers));
	if(temp > 1) % curve exceeds points, use clipped curve_t for reparametrization
		mask_left = []; curve_t_left = [];
	else % curve is shorter than points, use unclipped curve_t for reparametrization
		mask_left = idx(idx_curve{1}) == temp; % mask to unfiltered pts of first curve (not global curve_t)
		temp2 = idx_curve{1}(mask_left); % global index corrsponding to mask_left
		curve_t_left = t1(temp2) + t_offset_unclipped(temp2).*(t2(temp2)-t1(temp2));
	end
	temp = max(idx(inliers));
	if(temp < size(curve_pts,1)) % curve exceeds points, use clipped curve_t for reparametrization
		mask_right = []; curve_t_right = [];
	else % curve is shorter than points, use unclipped curve_t for reparametrization
		mask_right = idx(idx_curve{end}) == temp; % mask to unfiltered pts of first curve (not global curve_t)
		temp2 = idx_curve{end}(mask_right); % global index corrsponding to mask_right
		curve_t_right = t1(temp2) + t_offset_unclipped(temp2).*(t2(temp2)-t1(temp2));
	end

	% find extremal pts for each curve
	idx_boundary(:) = 0;
	if(num_curves == 1)
		t = curve_t(inliers); t(mask_left) = curve_t_left; t(mask_right) = curve_t_right; % t-parametrization of the curve with unrstricted ends (inliers only)
		t0 = min(t); t1 = max(t);
		curve_t(inliers) = (t-t0)./(t1-t0); % transforms t0->0 and t1->1
	else
		for i=1:num_curves
			if(~isempty(idx_curve{i}))
				% min
				if(i > 1) % corner
					[t, idx_min] = min(curve_t(idx_curve{i}));
					if(t > sampling_offsets(i,1)*1e-1) idx_boundary(i,1) = idx_curve{i}(idx_min); end % look only for pts 'inside' the curve, otherwise the curve length is not penalized (note: condition should be t>0 but soemtimes small ~1e-16 errors occur, note2: #bravo the relatively big offset 1e-1 (should be ~1e-10) is used to allow some oovershooting of the curve to delay the onset of the l2 penalty (empirically this helped with mini-oscillations, see #bravo))
				else % first curve - reparametrization of the left (free) end to match the endpts exactly in 1 iter
					t = curve_t(idx_curve{1}); t(mask_left) = curve_t_left; % t-parametrization of first curve (with possibly unclipped left end)
					t0 = min(t);
					if(abs(t0-1) > 1e-10)
						curve_t(idx_curve{1}) = (t-t0)/(1-t0); % transforms t0->0, 1->1
					end
				end

				% max
				if(i < num_curves)
					[t, idx_max] = max(curve_t(idx_curve{i}));
					if(t < 1-sampling_offsets(i,2)*1e-1) idx_boundary(i,2) = idx_curve{i}(idx_max); end % look only for pts 'inside' the curve, otherwise the curve length is not penalized (note: condition should be t<1 but soemtimes small ~1e-16 errors occur, note2: #bravo the relatively big offset 1e-1 (should be ~1e-10) is used to allow some oovershooting of the curve to delay the onset of the l2 penalty (empirically this helped with mini-oscillations, see #bravo))
				else % last curve - reparametrization of the right (free) end to match the endpts exactly in 1 iter
					t = curve_t(idx_curve{end}); t(mask_right) = curve_t_right; % t-parametrization of last curve (with possibly unclipped right end)
					t1 = max(t);
					if(abs(t1) > 1e-10)
						curve_t(idx_curve{end}) = t/t1; % transforms 0->0, t1->1
					end
				end
			else
				idx_boundary(i,:) = -1; % flag - indicates that corner loss applies for this curve, but cannot be bound to any actual input points belonging to this curve
			end
		end
	end

	% polish corner l2 loss - penalize only if both curves extend beyond their outermost points and bind to common points otherwise its acceptable situation and applying corner loss would slow down or stop correct convergence (also - binding the corner to two different points casuses problems - due to the condition t>0 etc for both points, the corner must be bound the the points that just corssed the threshold (during iterations) otherwise it's not stable and the whole curve starts to oscilate)
	for i=1:num_curves-1
		temp = i+[num_curves 1]; % linear indices in idx_boundary of the outermost points on the particular corner
		if(all(idx_boundary(temp))) % possibilities [x x], [x -1], [-1 -1]  where x > 0
			temp2 = find(idx_boundary(temp) > 0);
			if(numel(temp2) > 1)
				[~,temp3] = min(sum((control_pts(control_idx(i+1),:)-pts(idx_boundary(temp),:)).^2,2)); % distance from the corner pt (=first pt of i+1-th curve, which is simply the first control pt)
				idx_boundary(temp) = idx_boundary(temp(temp3));
			elseif(numel(temp2) == 1)
				idx_boundary(temp) = idx_boundary(temp(temp2));
			end % else both are -1, two curves without any points - no binding penalty...
		else % [x 0] or [0 0] - no binding penalty
			idx_boundary(temp) = 0;
		end
	end
	
% 	% DBG (pozn - nebude dobre fungovat zejmena kresleni pt matchingu - kvuli reparametrizaci to vypada ze projekce neni kolma)
% 	figure(1);
% 	plot(pts(:,2), pts(:,1), 'k+'); title(sprintf('iter=%d, diff=%.1e, loss=%.4e', iter, sqrt(diff2), sum(sqrt(dist2))));
% 	hold on;
% 	%plot(curve_pts(:,2), curve_pts(:,1), '-b'); % sampled pts
% 	for i=1:num_curves % plot corresponding curve
% 		pp = linspace(0,1,20).'.^(0:poly_order)*V*control_pts(control_idx(i):control_idx(i+1),:); % continuous curve
% 		plot(pp(:,2), pp(:,1), '-b');
% 		
% 		% asi bude fungovat blbe - nebylo opravene po pridani inliers
% 		%pp = (curve_t(idx_curve{i}).^(0:poly_order))*V*control_pts(control_idx(i):control_idx(i+1),:); % matching curve pt ofr each input pt
% 		%plot([pts(idx_curve{i}, 2).'; pp(:,2).'], [pts(idx_curve{i}, 1).'; pp(:,1).'], '-k');
% 		
% 		% sampling offsets
% 		%plot(curve_pts([curve_pts_idx(i) curve_pts_idx(i+1)-1],2), curve_pts([curve_pts_idx(i) curve_pts_idx(i+1)-1],1), 'om');
% 	end
% 	plot(pts(idx_boundary(idx_boundary>0),2), pts(idx_boundary(idx_boundary>0),1), 'om'); % boundary pts for each curve
% 	plot(pts(~inliers,2), pts(~inliers,1), 'or'); % outliers
% 	%plot([pts(:,2).'; curve_pts(idx,2).'], [pts(:,1).'; curve_pts(idx,1).'], '-r'); % closest sampled pt
% 	%plot([pts(:,2).'; curve_pts(idx+idx2,2).'], [pts(:,1).'; curve_pts(idx+idx2,1).'], '-g'); % neighboring sampled pt
% 	hold off;
% 	axis equal;
% 	drawnow;
% 	%pause;
	
	% construct linear system (gradient)
	w_irls = w./sqrt(dist2 + epsilon^2); % reweighting
	lhs(:) = 0; rhs(:) = 0; % reset system, build equations for each curve - these will overlap at shared control pts
	tv_corner = [V(1,:); sum(V,1)]; % corresponds to sampling at curve ends (t=0,1) - for corner loss that binds curve ends to input pts. These input pts are in idx_boundary
	cw = curvature_w(1)*sum(w(inliers));
	for i=1:num_curves
		span = control_idx(i):control_idx(i+1); % which indices in the linear system are affected by the processed curve
		m = idx_boundary(i,:) > 0; tv_temp = tv_corner(m, :); % mask of pts where corner-loss term applies
		tv = (curve_t(idx_curve{i}).^(0:poly_order))*V;
		lhs(span,span) = lhs(span,span) + tv.'*(tv.*w_irls(idx_curve{i})) + tv_temp.'*tv_temp + cw*curvature_term{i};
		rhs(span,:) = rhs(span,:) + tv.'*(pts(idx_curve{i}, :).*w_irls(idx_curve{i})) + tv_temp.'*pts(idx_boundary(i,m),:);
	end
	control_pts = lhs\rhs;
	
	% convergence - exit only if point matching is not required as output
	diff2 = max(sum((control_pts-control_old).^2,2)); % max between iteration displacement of any control point (squared)
	if(nargout == 1 && diff2 < abs_tol^2)
		break;
	end
end

% recalculate control pts to coeffs - return as cell
coeffs = cell(1,num_curves);
for i=1:num_curves
	coeffs{i} = V*control_pts(control_idx(i):control_idx(i+1), :);
end

% extra output args
if(nargout > 1)
	info = struct;
	info.control_pts = control_pts;
	info.curve_t = curve_t;
	info.curve_idx = curve_idx;
	info.dist2 = dist2;
	info.iter = iter;
	info.inliers = inliers;
	
	% some metrics
	if(poly_order > 1)
		% distance (arc-length) parametrization corresponding to each input pt projection
		curve_s = zeros(num_samples, 1);
		for i=1:num_curves
			curve_s(curve_pts_idx(i):curve_pts_idx(i+1)-1) = cumsum([boundary_pad; sqrt(sum((curve_pts(curve_pts_idx(i)+1:curve_pts_idx(i+1)-1,:)-curve_pts(curve_pts_idx(i):curve_pts_idx(i+1)-2,:)).^2,2))]); % approx curve distance at each sampled pt
		end
		info.curve_s = curve_s(idx+idx2).*t_offset + curve_s(idx).*(1-t_offset); % arc-length parametrization corresponding to each input pt (like curve_t but distance); linear interpolation of thw two nearest curve pts; note: may be slightly negative due to boundary-pad and other errors

		% length
		info.length = curve_s(curve_pts_idx(2:end)-1) + boundary_pad;
		% info.length = zeros(1,num_curves);
		% for i=1:num_curves
		% 	info.length(i) = sum(sqrt(sum((curve_pts(curve_pts_idx(i)+1:curve_pts_idx(i+1)-1,:)-curve_pts(curve_pts_idx(i):curve_pts_idx(i+1)-2,:)).^2,2)));
		% end

		% corner angle (dot product)
		info.angle = zeros(1,num_curves-1); % dot product (cos of angle) of unit derivatives at curve corner. 1=straight, -1=u-turn, 0=perpendicular
		if(num_curves > 1)
			dt1 = curve_dt(curve_pts_idx(1:end-1),:); dt1 = dt1./sqrt(sum(dt1.^2,2)); % unit derivative at curve beginning
			dt2 = curve_dt(curve_pts_idx(2:end)-1,:); dt2 = dt2./sqrt(sum(dt2.^2,2)); % unit derivative at curve end
			info.angle = sum(dt2(1:end-1,:).*dt1(2:end,:),2);
		end

		% max curvature
		info.max_curvature = zeros(1,num_curves);
		if(poly_order == 2) % max curvature for parabolas is easy to compute analytically
			for i=1:num_curves
				t = min(max(-(coeffs{i}(2,:)*coeffs{i}(3,:).')/(2*coeffs{i}(3,:)*coeffs{i}(3,:).'),0),1); % clipped parametrization of max curvature pt of the curve between [0,1] (only works at quadratic case)
				d1 = coeffs{i}(2,:)+2*t*coeffs{i}(3,:); % derivative at 't'
				info.max_curvature(i) = 2*abs(det([d1; coeffs{i}(3,:)]))/(d1*d1.')^(3/2);
			end
		else % numerical estimate for cubic and higher orders
			for i=1:num_curves
				% second derivatives at sampled pts
				d2t = (sampling_t(curve_pts_idx(i):curve_pts_idx(i+1)-1).^(0:poly_order-2).*(1:poly_order-1).*(2:poly_order))*coeffs{i}(3:end,:); % 2nd derivative at sampled pts
				dt = curve_dt(curve_pts_idx(i):curve_pts_idx(i+1)-1,:); % already sampled 1st derivative
				c = (dt(:,1).*d2t(:,2)-dt(:,2).*d2t(:,1))./sum(dt.^2,2).^(3/2); % signed curvature at sampled pts
				info.max_curvature2(i) = max(abs(c));
			end
		end
	else
		% arc-length parametrization and length
		info.length = sqrt(sum((control_pts(2:end,:)-control_pts(1:end-1,:)).^2,2));
		info.curve_s = curve_t.*info.length(curve_idx); % arc-length parametrization corresponding to each input pt (like curve_t but distance)
		
		% corner angle (dot product)
		info.angle = zeros(1,num_curves-1); % dot product (cos of angle) of unit derivatives at curve corner. 1=straight, -1=u-turn, 0=perpendicular
		if(num_curves > 1)
			dt = (control_pts(2:end,:)-control_pts(1:end-1,:))./info.length; % unit derivative of i-th curve
			info.angle = sum(dt(1:end-1,:).*dt(2:end,:),2);
		end

		% max curvature
		info.max_curvature = zeros(1,num_curves);
	end
end
end