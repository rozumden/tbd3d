
for i = 1:numel(seq)
	[latex_table,xx,yy] = estimate_speeds_tbd(frame{i}, r0_all(i));
	plot(xx,yy);
	[latex_table,xx,yy] = estimate_speeds(curves{i}, r0_all(i));
	plot(xx,yy);
end
