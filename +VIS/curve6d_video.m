function [] = curve6d_video(curves, szs, ind, matF, Vk, resz)
step_size = 0.5;
lw = 4;
tmin = curves(1).fit_iv(1);
tmax = curves(end).fit_iv(2) ;
modn = 10;

pnts = [];
for tind = tmin:step_size:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	pnt = [pnt interp1(ind,szs,tind)];
	if mod(tind,modn) == 0 && tind < tmax - 5
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
mltp = 10000;

%%% interpolate
pnts_use = [];
for tind = tmin:step_size:tmax
	pnt = evaluate_curve_coeff(curves, tind);
	pnt = ([pnt interp1(ind,szs,tind)]')./resz;
	pnt = mltp*from2Dto3D(pnt(1:2),pnt(3));
	kc = tind/modn;
	if floor(kc) < 1
    	draw_pnt = [pnt; pnts(4:7,1)];
    elseif ceil(kc) > size(pnts,2)
    	draw_pnt = [pnt; pnts(4:7,end)];
    else
    	p1 = pnts(4:7,floor(kc));
    	ckc = ceil(kc);
    	if ckc == kc
    		pdir = pnts(4:7,ckc);
	    	draw_pnt = [pnt; pdir];
    	else
	    	p2 = pnts(4:7,ckc);
	    	pdir = (1- (kc - floor(kc)))*p1 + (1 - (ckc - kc))*p2;
	    	draw_pnt = [pnt; pdir];
	    end
    end
    pnts_use = [pnts_use draw_pnt];
end
pnts_use(7,:) = pnts_use(7,:) / max(pnts_use(7,:));
% pnts_use(4:5,:) = pnts([5 4],:);
% pnts_use([2],:) = 3* pnts_use([2],:);
% pnts_use([3],:) = 5* pnts_use([3],:);
% pnts_use([1],:) = -3* pnts_use([1],:);
%%%%%
% pnts_use(1:3,:) = pnts_use([1 3 2],:);
% pnts_use(4:6,:) = pnts_use([5 6 4],:);
% pnts_use(4:5,:) = pnts_use([5 4],:);

last_pnt = [];
k = 1;
for tind = tmin:step_size:tmax
	pnt = pnts_use(:,k);
	k = k + 1;
	if isempty(last_pnt)
		plot3(pnt(1),pnt(2),pnt(3),'.g'); hold on
	else
		plot3([pnt(1) last_pnt(1)],[pnt(2) last_pnt(2)],[pnt(3) last_pnt(3)],'g','LineWidth',lw);
	end
	last_pnt = pnt;
end

% pbaspect([size(Vk,2)/size(Vk,1) 1 1]);
pbaspect([1 1 1]);
set(gca,'linewidth',3);
% xticks([0:50:size(Vk,1)]);
% yticks([0:50:size(Vk,2)]);

set(gca,'FontSize',15);
box on
ax = gca;
ax.BoxStyle = 'full';
% axis tight;
axis equal;
ylim([-30 70])
% axis normal
grid on;
view([-15,-85]);
% fig1 = gcf;

mlt_arr = 30;
for k = 1:size(pnts_use,2)
	if mod(k,20) ~= 1, continue; end
	sall = draw_line3(pnts_use(1:3,k),pnts_use(1:3,k)+pnts_use(4:6,k)*mlt_arr*pnts_use(7,k),...
                    'LineColor', [0.95 0.4 0.2],...
                    'LineWidth', 1 ,...
                    'ArrowDirection', 1,...
                    'ArrowLength', 3,....
                    'ArrowIntend', 2,...
                    'ArrowAngle', 45,...
                    'ArrowColor', [0.95 0.2 0.2]);
	alpha(sall,0); 
end

xticklabels([])
yticklabels([])
zticklabels([])

keyboard

writer = VID.VideoWriterWrapperDouble6D('~/projects/vis/', '6dplot', 'fps', 20);

gm = 0.5;
k = 0;
for tind = tmin:step_size:(tmax-50)
	k = k + 1;
	if mod(k,2) ~= 1, continue; end

	disp(tind);
    % set(0, 'CurrentFigure', fig2);
	frmi = ceil(tind/4);
	img = 255*(im2double(Vk(:,:,:,frmi)).^gm);
	zs = 255*ones(200, size(img,2), 3);
	img = cat(1, zs, img, zs);
    % imshow(IM); hold on
    % text(tx,4*ty,txt_orig,'FontSize',1.2*font_size,'Color',color_text);
    writer.write_buffer_img(img);
    % clf

    % set(0, 'CurrentFigure', fig1);
    draw_pnt = pnts_use(:,k);
	sall = draw_line3(draw_pnt(1:3),draw_pnt(1:3)+draw_pnt(4:6)*(mlt_arr-5)*draw_pnt(7),...
                    'LineColor', [0.95 0.4 0.2],...
                    'LineWidth', 1 ,...
                    'ArrowDirection', 1,...
                    'ArrowLength', 3,....
                    'ArrowIntend', 2,...
                    'ArrowAngle', 45,...
                    'ArrowColor', [0.95 0.2 0.2]);
	drawnow;
	inter = [100 250];
	alp = max(0,min(1, (tind-inter(1))/(inter(2) - inter(1))));
	% el = (1- alp)*(45) + alp*(-45);
	% view([el, 45]);
	el = (1- alp)*(-90) + alp*(-45);
	az = (1- alp)*(0) + alp*(-15);
	view([az,el]);

	grid on;
	writer.write();
	delete(sall);
end

writer.close();

keyboard