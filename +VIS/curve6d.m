function [] = curve6d(curves, szs, ind, matF)
step_size = 0.5;
lw = 4;
tmin = curves(1).fit_iv(1);
tmax = curves(end).fit_iv(2) + 1;

pnts = [];
for tind = tmin:step_size:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	pnt = [pnt interp1(ind,szs,tind)];
	if mod(tind,10) == 0 && tind < tmax - 5
		disp(['Done ' int2str(tind) ' \ ' int2str(tmax)]);
		itvr = find(ind > tind - 5 & ind < tind + 5);
		timestamps = ind(itvr);
		imgs = matF(:,:,:,itvr); 
		imgs_cell = mat2cell(imgs,size(imgs,1),size(imgs,2),size(imgs,3),repmat([1],1,numel(itvr))); imgs_cell = [imgs_cell(:)];
		[u v] = TD.findBallRotation(imgs_cell, timestamps);
		% usc = u*v;
		pnts = cat(1,pnts,[pnt u v]);
	end
end
pnts = pnts';
pnts_use = pnts;
hq = quiver3(pnts_use(1,:),pnts_use(2,:),pnts_use(3,:),pnts_use(4,:),pnts_use(5,:),pnts_use(6,:),0); 
hold on
hq.Color = [0.95 0.4 0.2];
hq.MaxHeadSize = 10;
hq.LineWidth = 2;


last_pnt = [];
for tind = tmin:step_size:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	pnt = [pnt interp1(ind,szs,tind)];
	if isempty(last_pnt)
		plot3(pnt(1),pnt(2),pnt(3),'.g'); hold on
	else
		plot3([pnt(1) last_pnt(1)],[pnt(2) last_pnt(2)],[pnt(3) last_pnt(3)],'g','LineWidth',lw);
	end
	last_pnt = pnt;
end

pnts_use(7,:) = pnts(7,:) / max(pnts(7,:));
pnts_use(6,:) = pnts(6,:);
for k = 1:size(pnts_use,2)
	draw_line3(pnts_use(1:3,k),pnts_use(1:3,k)+pnts_use(4:6,k)*30*pnts_use(7,k),...
                    'LineColor', [0.95 0.4 0.2],...
                    'LineWidth', 1 ,...
                    'ArrowDirection', 1,...
                    'ArrowLength', 3,....
                    'ArrowIntend', 2,...
                    'ArrowAngle', 45,...
                    'ArrowColor', [0.95 0.2 0.2]);
end