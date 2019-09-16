function [expf, ef] = estimate_exposure(frms)

ef = [];
for k = 1:numel(frms)-1
	frm = frms{k};
	frm1 = frms{k+1};
	if ~isempty(frm) && ~isempty(frm.Start) && ~isempty(frm1) && ~isempty(frm1.Start)
		f1 = norm(frm.Start - frm.End);
		f2 = norm(frm.Start - frm1.Start);
		ef = [ef f1/f2];
	end
end

ef = double(ef);
expf = mean(ef);
