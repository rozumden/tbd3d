function [ix] = frms2curve(fs, curve)
ix = zeros(2, numel(fs));

frms_all = [fs{:}];
with_direction = arrayfun(@(x) ~isempty(x.Direction), frms_all);
frms_all = frms_all(with_direction);

for frm = frms_all
	d1 = sum((curve - frm.Start').^2);
	[~, st_ind] = min(d1); st_ind = st_ind(1);

	d2 = sum((curve - frm.End').^2);
	[~, en_ind] = min(d2); en_ind = en_ind(1);
	ix(1, frm.instance) = st_ind;
	ix(2, frm.instance) = en_ind;
end

