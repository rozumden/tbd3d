fs = 35;
ns = [8 4 2 1];
load('../data/TbD-3D-n1_post.mat');
% clr1 = [0.85 0.47 0.32];
% clr2 = [0.1 0.25 0.89];
% clr3 = [0 0.4470 0.7410];
% clr4 = 'm';
% clr5 = [148 130 78]./255;

clr1 = [178, 56, 80]/255;
clr2 = [59, 139, 235]/255;
clr3 = [191, 227, 212]/255;
clr4 = [196, 219, 246]/255;
clr5 = [133, 144, 170]/255;

lw = 8;
bins = 240./ns;

barplot = true;

if false
	clf
	plot(bins, averages.tiou,'Color',clr1,'LineWidth',lw); hold on
	plot(bins, averages.tiou_nc,'Color',clr2,'LineWidth',lw);
	xlabel('fps');
	ylabel('TIoU')
	xlim([30 bins(end)]);
	xticks(bins);
	ylim([0.4 1])
	box off
	set(gca,'FontSize',fs)
	set(0,'defaulttextinterpreter','none')
	legend({'TbD','TbD-NC'});
	saveas(gcf,'~/tmp/tiou.png');
end

clf

if barplot
	b = bar([1 2 3 4], [averages.tiou3d; averages.tiou3d_nc; averages.nc3d3d; averages.tiou3d_nc3d_oracle]'); 
	b(1).FaceColor = clr4;
	b(2).FaceColor = clr5;
	b(3).FaceColor = clr1;
	b(4).FaceColor = clr2;
	xticklabels(bins);
else
	plot(bins, averages.tiou3d,'Color',clr4,'LineWidth',lw); 
	hold on
	plot(bins, averages.tiou3d_nc,'Color',clr5,'LineWidth',lw);
	plot(bins, averages.nc3d3d,'Color',clr1,'LineWidth',lw);
	plot(bins, averages.tiou3d_nc3d_oracle,'Color',clr2,'LineWidth',lw);
	xlim([30 bins(end)]);
	ylim([0.4 1])
	xticks(bins);
end

% plot(bins, averages.tiou3d_nc_oracle,'Color',clr4,'LineWidth',lw);
xlabel('fps');
ylabel('TIoU-3D')
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD','TbD-NC','TbD-3D','TbD-3D-O'},'Location','northwest');
saveas(gcf,'~/tmp/tiou3d.png');

if false
	clf
	plot(bins, averages.rerr,'b','LineWidth',lw); hold on
	plot(bins, averages.rerr_est,'c','LineWidth',lw);
	plot(bins, averages.rerr_gt,'m','LineWidth',lw);
	xlabel('fps');
	ylabel('Radius error')
	xlim([10 bins(end)]);
	box off
	set(gca,'FontSize',fs)
	set(0,'defaulttextinterpreter','none')
	legend({'TbD','TbD-3D','TbD-3D-O'});
	saveas(gcf,'~/tmp/errprs.png');
end
