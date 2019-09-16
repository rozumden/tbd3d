function [curve] = sequence_fit_part(frms, im)
Ts = zeros(size2(im));
frms_all = [frms{:}];
with_direction = arrayfun(@(x) ~isempty(x.Direction), frms_all);
frms_all = frms_all(with_direction);
	
for frm = frms_all
	if ~isempty(frm.bb)
		[H,~,~] = myTrajRender(size(frm.h), cellfun(@(x)fliplr(x).', frm.coeff, 'UniformOutput', false));
		% H = frm.h;
		[y,x] = find(H > 0);
		vals = H(sub2ind(size(H), y, x));
		inds = sub2ind(size(Ts), y+frm.bb(2)-1, x+frm.bb(1)-1);
		Ts(inds) = Ts(inds) + vals;
	else
		Ts = Ts + frm.h;
	end
end

[st, ~] = FIT.get_sten(frms_all(1));
[~, en] = FIT.get_sten(frms_all(end));
st = round(st); en = round(en);
% Tsf = imgaussfilt(Ts, 1);
[curve, hh] = FIT.fit_dp_seq(Ts, st, en, false);

