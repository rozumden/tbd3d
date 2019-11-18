% fs = 30;
ns = [8 4 2];
% load('../data/rotation_averages.mat');
load('../data/TbD-3D-rotations-abs.mat');

set(0,'defaulttextinterpreter','latex')

% lw = 10;
bins = 240./ns;
clf('reset')
hold on
if barplot
	if numel(ns) > 3
		b = bar([1 2 3 4], [rot_averages.err_u; rot_averages.err_or_u]'); 
		xticks([1 2 3 4]);
	else
		b = bar([1 2 3], [rot_averages.err_u(1:3); rot_averages.err_or_u(1:3)]'); 
		xticks([1 2 3]);
	end
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
ylim([0 45]);

h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', '~/tmp/err_u.pdf')

% saveas(gcf,'~/tmp/err_u.png');

clf('reset')

if barplot
	if numel(ns) > 3
		b = bar([1 2 3 4], rad2deg([rot_averages.err_rot; rot_averages.err_or_rot])'); 
	else
		b = bar([1 2 3], rad2deg([rot_averages.err_rot(1:3); rot_averages.err_or_rot(1:3)])'); 
	end
	b(1).FaceColor = clr1;
	b(2).FaceColor = clr2;
	xticklabels(bins);
else
	plot(bins, rot_averages.err_rot,'Color',clr1,'LineWidth',lw); hold on
	plot(bins, rot_averages.err_or_rot,'Color',clr2,'LineWidth',lw);
	xlim([30 bins(end)]);
	xticks(bins);
end

set(0,'defaulttextinterpreter','latex')

xlabel('fps');
ylabel('Angle Error [degrees]')
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});

h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', '~/tmp/err_rot.pdf')

saveas(gcf,'~/tmp/err_rot.png');
