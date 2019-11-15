function [szs_fixed] = depth_fit(inds, szs)
power_factor = 3;
expos = 1;
max_outliers = 1;

%% Split intervals
min_curve_length = 10;
wsize = 3;
szs_fixed = szs;

dc = [0; diff(szs)]; 
ddcind = zeros(size(szs));
for kki = 1:numel(inds)
	ind = inds(kki);
	if ind < wsize || ind > (inds(end)-wsize)
		continue;
	end
	nbr = szs(inds > ind-wsize & inds < ind+wsize);
	el = szs(kki);
	dfv = nbr - el;
	dfall = [0; diff(nbr)];
	split_pnt = kki - find(inds > ind-wsize,1) + 1;
	if numel(unique(sign(dfv))) == 3, continue; end %% not extreme value
	if max(abs(dfv)) <= 1.5, continue; end %% flat region
	%% perfectly flat on one side at least
	if max(abs(dfv(1:split_pnt-1))) < 1 || max(abs(dfv(split_pnt+1:end))) < 1, continue; end 
	s1 = unique(sign(dfall(1:split_pnt-1))); s1 = s1(s1 ~= 0);
	s2 = unique(sign(dfall(split_pnt+1:end))); s2 = s2(s2 ~= 0);
	if numel(s1) > 1 || numel(s2) > 1, continue; end %% left or right is not monotonous
	ddcind(kki) = 1;
end
%% remove small regions
for kki = find(ddcind)
	[~,last_ind] = find(ddcind(1:kki-1)); 
	if isempty(last_ind), continue; end
	last_ind = last_ind(end);
	if kki - last_ind < min_curve_length
		ddcind(last_ind:kki) = 1;
	end
end

%% remove clusters
dt = diff(bwdist(~ddcind'));
temp = [0 dt; dt 0];
break_points = ismember(temp',[1 -1], 'rows')' | ismember(temp',[1 0], 'rows')';

break_points(1) = 1;
break_points(end) = 1;
break_points = find(break_points);
for kki = 2:numel(break_points)
	
	st_fit = inds(break_points(kki-1));
	en_fit = inds(break_points(kki));

	max_power = max(2, min(6, round((en_fit-st_fit)/power_factor) ));
	if en_fit-st_fit == 0, max_power = 1; end
	inds_inl = inds(break_points(kki-1):break_points(kki));
	szs_inl = szs(break_points(kki-1):break_points(kki));
	coeff = depth_lsq_fit(inds_inl, szs_inl, max_power, expos);
	for k2 = break_points(kki-1):break_points(kki)
		szs_fixed(k2) = evaluate_coeff(coeff{1}, inds(k2));
	end
end



function [coeff] = depth_lsq_fit(inds, szs, max_power, expos)
A = [];
b = [];
for ifit = 1:numel(inds)
	Ast = [];
	for pow1 = 0:max_power
		Ast=[Ast [inds(ifit).^pow1;]];
 	end
 	
	A = [A; Ast];
	b = [b; szs(ifit)'];

end

if isempty(A)
	coeff = {ones(2,2)};
	return;
end

xs = A \ b;
nvars = (max_power+1);
options = optimset('Display','none');
xs = lsqlin(A,b,[],[],[],[],-Inf*ones(nvars,1),Inf*ones(nvars,1),[],options);

coeff = {reshape(xs, [1 max_power+1])};


