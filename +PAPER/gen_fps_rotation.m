% fs = 30;
ns = [8 4 2 1];
% load('../data/rotation_averages.mat');
load('../data/TbD-3D-rotations-abs.mat');

% lw = 10;
bins = 240./ns;
clf('reset')
hold on
if barplot
	b = bar([1 2 3 4], [rot_averages.err_u; rot_averages.err_or_u]'); 
	b(1).FaceColor = clr1;
	b(2).FaceColor = clr2;
	xticklabels(bins);
else
	plot(bins, rot_averages.err_u,'Color',clr1,'LineWidth',lw);
	plot(bins, rot_averages.err_or_u,'Color',clr2,'LineWidth',lw);
	xlim([30 bins(end)]);
	xticks(bins);
end
xlabel('fps');
ylabel('Axis Error [degrees]')
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_u.png');
clf('reset')

if barplot
	b = bar([1 2 3 4], rad2deg([rot_averages.err_rot; rot_averages.err_or_rot])'); 
	b(1).FaceColor = clr1;
	b(2).FaceColor = clr2;
	xticklabels(bins);
else
	plot(bins, rot_averages.err_rot,'Color',clr1,'LineWidth',lw); hold on
	plot(bins, rot_averages.err_or_rot,'Color',clr2,'LineWidth',lw);
	xlim([30 bins(end)]);
	xticks(bins);
end

xlabel('fps');
ylabel('Angle Error [degrees]')
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_rot.png');
