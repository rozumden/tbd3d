
est_rot = repmat({[]}, 1, numel(seq));

% for i = 1:numel(seq)
% for i = seq_id
parfor i = 1:numel(seq)
	disp(seq(i).name);
	warning('off','all');
	
	% imgs = num2cell(matF_hs_cell{i},[1 2 3]); 
	% imgs = [imgs(:)];
	imgs = matF_hs_cell{i};

	matF_est_WB = matF_cell{i};
	for kkk = 1:numel(matF_est_WB)
		for k = 1:3, matF_est_WB{kkk}(:,:,k) = matF_est_WB{kkk}(:,:,k) ./ fc(k); end
		matF_est_WB{kkk} = ((matF_est_WB{kkk}.*reshape(WB,1,1,[])/(max(WB))).^gamma_coef);
	end

	matF_gt_WB = matF_gt_cell{i};
	for kkk = 1:numel(matF_gt_WB)
		for k = 1:3, matF_gt_WB{kkk}(:,:,k) = matF_gt_WB{kkk}(:,:,k) ./ fc(k); end
		matF_gt_WB{kkk} = ((matF_gt_WB{kkk}.*reshape(WB,1,1,[])/(max(WB))).^gamma_coef);
	end

	batch = 32;
	ki = 0;
	est_rot{i} = struct('hs',[],'or',[],'est',[]);

	% eval_at = [batch:batch:numel(imgs)];
	eval_at = [batch:4:1.5*batch];

	for k = eval_at
		disp([int2str(ki) ' : ' seq(i).name]);
		ki = ki + 1;

		iv = k-batch+1:k;
		imgs_use = imgs(iv);
		res_hs = [];
		[res_hs.u res_hs.rot] = TD.findBallRotation(imgs_use);
		est_rot{i}(ki).hs = res_hs;

		in_resl = (ind_gt{i}-1)*n+1;
		inl_gt = in_resl >= iv(1) & in_resl < iv(end)+1;
		ind_gt_use = in_resl(inl_gt);
		% matF_gt_use = num2cell(matF_gt_WB(:,:,:,inl_gt),[1 2 3]); matF_gt_use = [matF_gt_use(:)];
		matF_gt_use = matF_gt_WB(inl_gt);
		res_gt = [];
		[res_gt.u res_gt.rot] = TD.findBallRotation(matF_gt_use,ind_gt_use);
		est_rot{i}(ki).or = res_gt;
		
		in_resl = (ind_est{i}-1)*n+1;
		inl_est = in_resl >= iv(1) & in_resl < iv(end)+1;
		ind_est_use = in_resl(inl_est);
		% matF_use = num2cell(matF_est_WB(:,:,:,inl_est),[1 2 3]); matF_use = [matF_use(:)];
		matF_use = matF_est_WB(inl_est);
		res = [];
		[res.u res.rot] = TD.findBallRotation(matF_use, ind_est_use);
		est_rot{i}(ki).est = res;
		
	end

end

err_u = zeros(1,numel(seq));
err_rot = zeros(1,numel(seq));

err_or_u = zeros(1,numel(seq));
err_or_rot = zeros(1,numel(seq));



std_u = zeros(1,numel(seq));
std_rot = zeros(1,numel(seq));

std_or_u = zeros(1,numel(seq));
std_or_rot = zeros(1,numel(seq));


u_y = zeros(2,0);
rot_y = zeros(2,0);
for i = 1:numel(seq)
	nel = numel(est_rot{i});
	ierr_u = zeros(1,nel);
	ierr_or_u = zeros(1,nel);

	rot_hs = zeros(1,nel);
	rot_or = zeros(1,nel);
	rot_est = zeros(1,nel);

	for ki = 1:nel
		ierr_u(ki) = est_rot{i}(ki).est.u(:).'*est_rot{i}(ki).hs.u(:);
		ierr_or_u(ki) = est_rot{i}(ki).or.u(:).'*est_rot{i}(ki).hs.u(:);

		rot_hs(ki) = est_rot{i}(ki).hs.rot;
		rot_or(ki) = est_rot{i}(ki).or.rot;
		rot_est(ki) = est_rot{i}(ki).est.rot;
	end
	err_u(i) = mean(rad2deg(acos(abs(ierr_u))));
	err_or_u(i) = mean(rad2deg(acos(abs(ierr_or_u))));
	err_rot(i) = mean(abs(rot_hs - rot_est));
	err_or_rot(i) = mean(abs(rot_hs - rot_or)); 

	std_u(i) = std(rad2deg(acos(abs(ierr_u))));
	std_or_u(i) = std(rad2deg(acos(abs(ierr_or_u))));
	std_rot(i) = std(abs(rot_hs - rot_est));
	std_or_rot(i) = std(abs(rot_hs - rot_or));

	u_y = [u_y [rad2deg(acos(abs(ierr_u))); rad2deg(acos(abs(ierr_or_u)))]];
	rot_y = [rot_y [abs(rot_hs - rot_est); abs(rot_hs - rot_or)]];
end

u_p = anova1(u_y');
rot_p = anova1(rot_y');

if true
	fprintf('[Oracle] Mean angle of axis error is %.4f\n', mean(err_or_u));
	fprintf('[Oracle] Std angle of axis error is %.4f\n', mean(std_or_u));
	fprintf('Mean angle of axis error is %.4f\n', mean(err_u));
	fprintf('Std angle of axis error is %.4f\n', mean(std_u));
	fprintf('[Oracle] Mean absolute error of angular speed is %.4f\n', mean(err_or_rot));
	fprintf('[Oracle] Std absolute error of angular speed is %.4f\n', mean(std_or_rot));
	fprintf('Mean absolute error of angular speed is %.4f\n', mean(err_rot));
	fprintf('Std absolute error of angular speed is %.4f\n', mean(std_rot));
end

