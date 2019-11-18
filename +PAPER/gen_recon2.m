
is = [2 5 5 4];
ks = [1 31 63 3];
n = 8;

sz = 5;
for kkk = 1:numel(is)
	i = is(kkk);
	k = floor(ind_gt{i}(ks(kkk)));

	scl = 4;
	new_coeff = gt_coeffs{i}{k};
	for ci = 1:n, new_coeff{ci} = scl*new_coeff{ci}; end
	% new_coeff = [gt_coeffs{i}{1}{1}(:,1) gt_coeffs{i}{1}{end}(:,1)];
	% new_coeff = {[new_coeff(:,1) new_coeff(:,2)-new_coeff(:,1)]*scl};
	inp = imresize(Vk_WB{i}(:,:,:,k),scl);
	if kkk == 1
		new_coeff = new_coeff(1:n-1);
	end
	Hall = myTrajRender(size2(inp), new_coeff, [0 1]);

	clr = [0 255 0];
	for ui = 1:3
		msk = logical(zeros(size(inp)));
		msk(:,:,ui) = conv2(Hall,ones(sz),'same') > 0;
		inp(msk) = clr(ui); 
	end

	for tk = 4
		clr = [255 0 0];
		position = round(new_coeff{tk}(:,1) + 0.5*new_coeff{tk}(:,2));
		for ii = -sz:sz
			for jj = -sz:sz
				for ui = 1:3, inp(position(2)+ii,position(1)+jj,ui) = clr(ui); end
			end
		end
	end
	HM = conv2(Hall,ones(scl*size(matM_hs{i}(:,:,:,1))),'same');
	[yy,xx]=find(HM>0);
	xs = [min(xx(:)) max(xx(:))];
	df = xs(2) - xs(1);
	th = 450;
	if df < th
		dfth = th - df;
		xs(1) = xs(1) - dfth/2;
		xs(2) = xs(2) + dfth/2;
	end

	ys = [min(yy(:)) max(yy(:))];
	df = ys(2) - ys(1);
	th = 400;
	if df < th
		dfth = th - df;
		ys(1) = ys(1) - dfth/2;
		ys(2) = ys(2) + dfth/2;
	end

	inp = inp(ys(1):ys(2),xs(1):xs(2),:);
	imwrite(inp,['~/img_' int2str(kkk) 'use.png']);
end