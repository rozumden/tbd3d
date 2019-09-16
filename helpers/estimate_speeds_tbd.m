function [srctable,xx,yy] = estimate_speeds_tbd(frms, r0)
srctable = '';
xx = []; yy = [];
for k = 1:numel(frms)
	frm = frms{k};
	if isempty(frm)
		continue;
	end
	
	for ti = 0:0.25:0.99
		spd = norm(evaluate_dcoeff(frm.coeff{1},ti)) / r0;

		srctable = [srctable sprintf('(%.2f, %.4f) ',k+ti, spd)];
		xx = [xx k+ti]; yy = [yy spd];
	end
end

fid = fopen('~/projects/vis/report/report.tex','wt');
fprintf(fid, srctable);
fclose(fid);