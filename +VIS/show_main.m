function [] = show_main(im0,frame,n)
% figure(1);
clf;
image(im0);
hold on;
lw = 10;
for kk = 0:min(0,n-1)
	for kkk = 1:numel(frame{n-kk})
		traj = frame{n-kk}(kkk).TrajectoryXY;
		if strcmp(frame{n-kk}(kkk).caseused, 'FMOd')
			% frame{n-kk}(kkk).TrajectoryXY = traj(:,1:2:end);
		elseif strcmp(frame{n-kk}(kkk).caseused, 'Std')
			rectangle('Position',frame{n-kk}(kkk).bb);
		end
		if kk == 0,
			% 1/((kk+1)^2)
			frame{n-kk}(kkk).show;
		else
			frame{n-kk}(kkk).show_traj;
		end
		frame{n-kk}(kkk).TrajectoryXY = traj;
	end
end
drawnow;
