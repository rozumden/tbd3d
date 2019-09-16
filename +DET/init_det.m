function frm = init_det(detec, IM)
frm = detec.detect(IM);

if ~isempty(frm) && ~isempty(frm(1).fittingScore)
	[~,ind] = max([frm.fittingScore]);
	% [~,ind] = max([frm.Area]);
	if isempty(ind)
		return
	end
	frm = frm(ind);
end
