
ns = [8 4 2 1];
load('../data/rotation_averages.mat');

lw = 10;
bins = 240./ns;
figure;
hold on
plot(bins, rot_averages.err_u,'b','LineWidth',lw);
plot(bins, rot_averages.err_or_u,'g','LineWidth',lw);
xlabel('fps');
ylabel('Average Angle')
xlim([10 bins(end)]);
box off
set(gca,'FontSize',20)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_u.png');
clf

plot(bins, rot_averages.err_rot,'b','LineWidth',lw);
plot(bins, rot_averages.err_or_rot,'g','LineWidth',lw);
xlabel('fps');
ylabel('Average Error [%]')
xlim([10 bins(end)]);
box off
set(gca,'FontSize',20)
set(0,'defaulttextinterpreter','none')
legend({'TbD-3D','TbD-3D-O'});
saveas(gcf,'~/tmp/err_rot.png');
clf