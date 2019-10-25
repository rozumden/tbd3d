fs = 30;
ns = [8 4 2 1];
load('../data/TbD-3D-n1_post.mat');

lw = 10;
bins = 240./ns;
clf
plot(bins, averages.tiou,'b','LineWidth',lw); hold on
plot(bins, averages.tiou_nc,'g','LineWidth',lw);
xlabel('fps');
ylabel('TIoU')
xlim([10 bins(end)]);
ylim([0.4 1])
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD','TbD-NC'});
saveas(gcf,'~/tmp/tiou.png');
clf

plot(bins, averages.tiou3d,'b','LineWidth',lw); hold on
plot(bins, averages.tiou3d_nc,'g','LineWidth',lw);
plot(bins, averages.nc3d3d,'c','LineWidth',lw);
plot(bins, averages.tiou3d_nc_oracle,'r','LineWidth',lw);
plot(bins, averages.tiou3d_nc3d_oracle,'m','LineWidth',lw);
xlabel('fps');
ylabel('TIoU-3D')
xlim([10 bins(end)]);
ylim([0.4 1])
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD','TbD-NC','TbD-3D','TbD-O','TbD-3D-O'});
saveas(gcf,'~/tmp/tiou3d.png');
clf

if false
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
