function [curve, H] = fit_dp_seq(H, st, en, vis)
if nargin < 3
	vis = false;
end
th = 0.1;
sigma_hist = 0;
step_size = 2;
sigma_psf = 2;

[curve, vals] = fit_dp_start_seq(H, st, en, th, sigma_hist, step_size, sigma_psf);
[curve_t, vals_t] = fit_dp_start_seq(H', st([2 1]), en([2 1]), th, sigma_hist, step_size, sigma_psf);	

if size(curve_t,2) > size(curve,2)
	curve = curve_t([2 1],:);
	vals = vals_t;
end

if vis
	Hvis = FIT.make_image(H / max(H(:)), curve, curve);
	figure; imshow(Hvis);
end

function [curve, vals] = fit_dp_start_seq(H, st, en, th, sigma_hist, step_size, sigma_psf)
if ~exist('sigma_hist')
	sigma_hist = 0;
end
if ~exist('th')
	th = 0.1;
end
if ~exist('step_size')
	step_size = 1;
end
[w1,h1] = size(H);
turn_eps = 0;
turn_bias = 0.1;
start_trun = 10;
flipped = false;

if en(1) < st(1) 
	st(1) = h1 - st(1) + 1;
	en(1) = h1 - en(1) + 1;
	H = fliplr(H);
	flipped = true;
end 

Hest = H; Hest(Hest < 0) = 0;
if sigma_psf > 0, Hest = imgaussfilt(Hest, sigma_psf); end	

% H = H / sum(sum(H(H>0)));
H = H / max(H(:));
J = zeros(w1,h1);
P = zeros(w1,h1);

% start_eps = quantile([Hest(Hest>0)],0.9);
% end_eps = quantile([Hest(:)],0.8);
% end_eps = max(Hest(:))*0.22;
if false
	start_eps = 0.1;
	end_eps = 0.1;

	% prior probability of start
	S = ([1:h1]-st(1));
	S = sign(S) .* (S.^2);
	S(S > start_trun) = start_trun;
	S = S*start_eps;

	% prior probability of end
	E = (en(1)-[1:h1]);
	% E = sign(E) .* (E.^2);
	% E(E > start_trun) = start_trun;
	E = repmat(E*end_eps, [w1 1]);
else
	S = -Inf*ones(w1,h1);
	S(round(st(2)),round(st(1))) = 0;
	E = -Inf*ones(w1,h1);
	E(round(en(2)),round(en(1))) = 0;
end

J(:,1) = H(:,1)+S(:,1);
P(:,1) = [1:w1];

nsz = 2*step_size + 1;
for i = 2:h1
	%% add data term
	if step_size == 1
		mat = [ [-Inf; J(1:end-1,i-1)] J(:,i-1) [J(2:end,i-1); -Inf]]';
	elseif step_size == 2
		mat = [ [-Inf; -Inf; J(1:end-2,i-1)] ...
			  [-Inf; J(1:end-1,i-1)] J(:,i-1) [J(2:end,i-1); -Inf] ...
			  [J(3:end,i-1); -Inf; -Inf]]';
	else
		error('Unknown step size');
	end
	%% add start penalty term
	mat = [mat; [repmat(S(:,i),[1 w1])]];

	%% add direction penalty term
	direc = P(:,i-1) - [1:w1]';
	direc(1) = min([0 direc(1)]);
	direc(2) = min([1 direc(2)]);
	direc(end) = max([0 direc(end)]);
	direc(end-1) = max([-1 direc(end-1)]);
	inds = [1:w1]' - direc;
	inds = inds(P(:,i-1) > 0); direc = direc(P(:,i-1) > 0);
	direc_mat = repmat((turn_bias*(-([step_size:-1:0 1:step_size].^2)'-1)), [1 w1]);
	direc_mat(sub2ind(size(direc_mat), direc+step_size+1, inds)) = turn_eps;
	mat(1:nsz,:) = mat(1:nsz,:) + direc_mat;

	[J(:,i), ind] = max(mat);
	J(:,i) = J(:,i) + H(:,i);
	ind = ind - step_size - 1;
	P(:,i) = [1:w1]+ind;
	P(ind > step_size, i) = 0; 
end

%% add end penalty term
J = J + E;

[~, ii] = find(J == max(J(:))); ii = min(ii);

traj = zeros(1, ii);
[~, traj(end)] = max(J(:,ii));
for i = ii-1:-1:1
	if traj(i+1) == 0, break; end %% start
	traj(i) = P(traj(i+1), i+1);
end

xs = [1:ii];
curve = [xs(traj>0); traj(traj>0)];

vals = H(sub2ind(size(H),curve(2,:),curve(1,:)));

if flipped
	% vals = fliplr(vals);
	% curve = fliplr(curve);
	curve(1,:) = h1 - curve(1,:) + 1;
end

