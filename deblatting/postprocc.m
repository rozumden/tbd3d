function [coeff, len, pnts] = postprocc(coeff, bb, parts, ivl, refl)
if nargin < 4
	ivl = [0 1];
end
if nargin < 5
	refl = true;
end
for c = 1:numel(coeff)
	if refl
		coeff{c} = fliplr(coeff{c}).';
	end
	if numel(bb) == 2
		coeff{c}(:,1) = coeff{c}(:,1) + bb;
	end
end
% dof = 6;
% w = cellfun(@(x)size(x,2),coeff) < dof;
% coeff(w) = cellfun(@(x)[x zeros(2,dof-size(x,2))], coeff(w), 'UniformOutput', false);


len = zeros(1,numel(coeff));
for c = 1:numel(coeff)
	max_power = (numel(coeff{c}) / 2) - 1;
	ti = linspace(ivl(1),ivl(2),parts); % linearize the curve by this many segments
	t = [ones(length(ti),1)];
	for pow1 = 1:max_power
		t = [t ti(:).^pow1];
	end
    p = t*coeff{c}.'; %'
    len(c) = sum(sqrt(sum((p(2:end,:)-p(1:end-1,:)).^2,2)));
end
fracs = len / sum(len);
fracs(isnan(fracs)) = 0;

parts0 = floor(parts * fracs);
parts0(end) = parts - sum(parts0(1:end-1));
pnts = [];

for c = 1:numel(coeff) 
	max_power = (numel(coeff{c}) / 2) - 1;
	dt = linspace(ivl(1),ivl(2),parts0(c));
	ts = [ones(length(dt),1)];
	for pow1 = 1:max_power
		ts = [ts dt(:).^pow1];
	end
	pnts = [pnts; ts*coeff{c}.'];
end
if numel(bb) == 8
	pnts = TRACK.change_coor(bb, pnts(:,2), pnts(:,1))';
end