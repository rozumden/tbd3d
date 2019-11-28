function [] = generate_imgs(matF, matM, frame, matF_hs, matM_hs, ind_gt, iv, names)
if ~exist('names','var')
	names = [];
end
% if nargin < 4
% 	iv = 41:0.5:43;
% end
fc = [1.9 1 1.8];
WB = [2 1 2]; gamma_coef = 0.4;
for k = 1:3
	matF(:,:,k,:) = matF(:,:,k,:) ./ fc(k);
end

matF = ((matF.*reshape(WB,1,1,[])/(max(WB))).^gamma_coef);

th = 0.7;
draw_rotation = true;  %% draw rotation axis

n = size(matF_hs,4) / numel(frame);
frmt = '.png';

if draw_rotation
	% itvr = [(iv-1)*n + 1];
	itvr = 1:40;
	r = ceil(sqrt(sum(sum(matM_hs(:,:,:,1)))/pi))-1;

	imgs = matF(:,:,:,itvr); imgs = imgs(end/2-r+1:end/2+r,end/2-r+1:end/2+r,:,:);
	imgs_cell = mat2cell(imgs,size(imgs,1),size(imgs,2),size(imgs,3),repmat([1],1,numel(itvr))); imgs_cell = [imgs_cell(:)];
	timestamps = ind_gt(itvr);
	[u v] = TD.findBallRotation(imgs_cell, timestamps,1.5);
	
	imgs_gt = matF_hs(:,:,:,itvr); imgs_gt = imgs_gt(end/2-r+1:end/2+r,end/2-r+1:end/2+r,:,:);
	imgs_gt_cell = mat2cell(imgs_gt,size(imgs_gt,1),size(imgs_gt,2),size(imgs_gt,3),repmat([1],1,numel(itvr))); imgs_gt_cell = [imgs_gt_cell(:)];
	timestamps = ind_gt(itvr);
	[ugt vgt] = TD.findBallRotation(imgs_gt_cell, timestamps);

	% plot(2*r*(i-1) + r+sign(ugt(3))*ugt(2)*r-1, 2*r*(j-1)+ r+sign(ugt(3))*ugt(1)*r-1,'.g')
	% plot(2*r*(i-1) + r+sign(u(3))*u(2)*r-1, 2*r*(j-1)+ r+sign(u(3))*u(1)*r-1,'.r')
end

ki = 0;
for i = iv
	ki = ki + 1;
	if false
		it = (i-1)*n + 1;
	else
		[~,it] = find(ind_gt == i);
	end
	F3D = matF(:,:,:,it); 
	MM = matM(:,:,:,it); 

	GT = matF_hs(:,:,:,(i-1)*n + 1); 
	GTM = matM_hs(:,:,:,(i-1)*n + 1); 
	
	if true %% make smaller
		r = ceil(sqrt(sum(GTM(:))/pi))-1;
		GTM = GTM(end/2-r+1:end/2+r,end/2-r+1:end/2+r);
		GT = GT(end/2-r+1:end/2+r,end/2-r+1:end/2+r,:);
		F3D = F3D(end/2-r+1:end/2+r,end/2-r+1:end/2+r,:).*GTM;
		MM = MM(end/2-r+1:end/2+r,end/2-r+1:end/2+r,:).*GTM;
	end

	if draw_rotation
		i1 = repmat(GTM < 1-th,[1 1 3]); 
		F3D(i1) = 1-F3D(i1);
		GT(i1) = 1-GT(i1);

		if false
			clr = [255 50 0];
			for ui = 1:3
				for ii = -2:2
					for jj = -2:2
						if ii*jj ~= 0, continue; end
						GT(round(r+sign(ugt(3))*ugt(1)*r-1)+ii,round(r+sign(ugt(3))*ugt(2)*r-1)+jj,ui) = clr(ui);
						F3D(round(r+sign(u(3))*u(1)*r-1)+ii,round(r+sign(u(3))*u(2)*r-1)+jj,ui) = clr(ui);
					end
				end
			end 
		end

		imwrite(GT,['~/tmp/recon/gt_' int2str(it) frmt]);
		% imwrite(GTM,['~/gtM_' int2str(it) frmt]);
		imwrite(F3D,['~/tmp/recon/tbd3d_' int2str(it) frmt]);
		% imwrite(MM,['~/tbd3dM_' int2str(it) frmt]);
		% if round(i) == i
		% 	imwrite(frame{i}.f.^(1/2.2),['~/tbd_' int2str(it) frmt]);
		% 	imwrite(mean(frame{i}.f,3),['~/tbdM_' int2str(it) frmt]);
		% end
	else
		fr = frame{round(i)};
		% FMO = fr.f.^(1/2.2);
		% FMOM = mean(FMO.^(2.2),3);
		[f_img, FMOM] = estimateFM_motion_template(fr.im_c, fr.bgr_c, fr.h, [], zeros(size2(fr.f)), [], []);
		FMO = f_img.^(1/2.2);
		if size(FMO,1) < size(GT,1)
			psize = (size(GT,1) - size(FMO,1))/2;
			hlf = psize == round(psize);
			psize = floor(psize);
			FMO = padarray(FMO,[psize psize]);
			FMOM = padarray(FMOM,[psize psize]);
			if ~hlf
				FMO = padarray(FMO,[1 1],'post');
				FMOM = padarray(FMOM,[1 1],'post');
			end
		elseif  size(FMO,1) > size(GT,1)
			xCenter = round(size(FMO,1)/2);
			radius = floor(size(GT,2)/2);
			xLeft = xCenter - radius;
			FMO = imcrop(real(FMO), [xLeft, xLeft, size(GT,1)-1, size(GT,2)-1]);
			FMOM = imcrop(real(FMOM), [xLeft, xLeft, size(GT,1)-1, size(GT,2)-1]);
		end

		if size(F3D,1) < size(GT,1)
			psize = (size(GT,1) - size(F3D,1))/2;
			hlf = psize == round(psize);
			psize = floor(psize);
			F3D = padarray(F3D,[psize psize]);
			MM = padarray(MM,[psize psize]);
			if ~hlf
				F3D = padarray(F3D,[1 1],'post');
				MM = padarray(MM,[1 1],'post');
			end
		elseif size(F3D,1) > size(GT,1)
			xCenter = round(size(F3D,1)/2);
			radius = floor(size(GT,2)/2);
			xLeft = xCenter - radius;
			F3D = imcrop(real(F3D), [xLeft, xLeft, size(GT,1)-1, size(GT,2)-1]);
			MM = imcrop(real(MM), [xLeft, xLeft, size(GT,1)-1, size(GT,2)-1]);
		end

		if true
			MM = 1-MM; GTM = 1-GTM; FMOM = 1-FMOM;
			i1 = repmat(MM > th,[1 1 3]); F3D(i1) = 1-F3D(i1);
			i2 = repmat(GTM > th,[1 1 3]); GT(i2) = 1-GT(i2);
			i3 = repmat(FMOM > th,[1 1 3]); FMO(i3) = 1-FMO(i3);
		end
		hmask = conv2(fr.Ti,imresize(GTM,1.4),'same') > 0;
		[y,x] = find(hmask);
		x = sort(x(:)); y = sort(y(:));
		inimg = imresize(fr.im_c(y(1):y(end),x(1):x(end),:).^(1/2.2), [size(GTM,1) NaN]);

		% inimg = imresize(fr.im_c(20:end-10,20:end-20,:).^(1/2.2), [size(GTM,1) NaN]);

		% jnt = cat(2,FMO,FMOM(:,:,[1 1 1]),F3D,MM(:,:,[1 1 1]),GT,GTM(:,:,[1 1 1]));
		jnt = cat(2,inimg,F3D,MM(:,:,[1 1 1]),GT,GTM(:,:,[1 1 1]));
		% jnt = imresize(jnt, 100/size(inimg,2));
		if isempty(names)
			imwrite(jnt,['~/joint_' int2str(it) frmt]);
		else
			% imwrite(jnt,['~/joint_' names{ki} frmt]);
			imwrite(inimg,['~/img_' names{ki} frmt]);
			% imwrite(F3D,['~/est_' names{ki} frmt]);
			% imwrite(MM,['~/estM_' names{ki} frmt]);
			% imwrite(GT,['~/gt_' names{ki} frmt]);
			% imwrite(GTM,['~/gtM_' names{ki} frmt]);
		end			
	end
end
