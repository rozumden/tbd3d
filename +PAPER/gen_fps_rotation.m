fs = 30;
ns = [8 4 2 1];
load('../data/rotation_averages.mat');

lw = 10;
bins = 240./ns;
clf
hold on
plot(bins, rot_averages.err_u,'b','LineWidth',lw);
plot(bins, rot_averages.err_or_u,'g','LineWidth',lw);
xlabel('fps');
ylabel('Average angle bw. axes [degrees]')
xlim([10 bins(end)]);
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_u.png');
clf

plot(bins, 100*rot_averages.err_rot,'b','LineWidth',lw); hold on
plot(bins, 100*rot_averages.err_or_rot,'g','LineWidth',lw);
xlabel('fps');
ylabel('Average error [%]')
xlim([10 bins(end)]);
box off
set(gca,'FontSize',fs)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_rot.png');
