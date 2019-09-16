function [curve, H] = fit_dp(H, st, en, vis)
if nargin < 3
	vis = false;
end
th = 0.1;
sigma_hist = 0;
step_size = 2;
sigma_psf = 2;

% if sigma_psf > 0, H = imgaussfilt(H, sigma_psf); end	


[curve, vals, J, P] = FIT.fit_dp_start(H, st, en, th, sigma_hist, step_size, sigma_psf);
[curve_t, vals_t, J_t, P_t] = FIT.fit_dp_start(H', st([2 1]), en([2 1]), th, sigma_hist, step_size, sigma_psf);	

if size(curve_t,2) > size(curve,2)
	curve = curve_t([2 1],:);
	vals = vals_t;
	J = J_t';
	P = P_t';
end

mrsz = 20;
Z = logical(zeros(size(P)));
if vis
	% figure; imshow(permute(fr.im_c,[2 1 3])); hold on; plot(a1:a2,traj(a1:a2),'r'); 
	% figure; imshow(Hv); hold on; plot(1:h, traj, 'b'); plot(a1:a2,traj(a1:a2),'r'); 
	% figure; plot(1:numel(vals), vals);
	Hvis = FIT.make_image(H / max(H(:)), curve, curve);
	imshow(Hvis); hold on;
	plot(st(1),st(2),'xg','MarkerSize',mrsz,'LineWidth',mrsz);
	plot(en(1),en(2),'xy','MarkerSize',mrsz,'LineWidth',mrsz);
	saveas(gca, 'dp_output.png');
	% imwrite(Hvis, 'dp_output.png','Mode','lossless');
	J = J - min(J(:)); J = 1 - J/max(J(:));
	% J(P == 0) = 1;
	% Jvis = FIT.make_image(J, curve, curve, st, en);
	Jvis = J;
	imshow(Jvis); 
	imwrite(Jvis, 'dp_energy.png','Mode','lossless');
	P(:,1) = 0;
	% Pvis = FIT.make_image(double(P == 0), curve, curve, st, en);
	Pvis = repmat(double(P == 0), [1 1 3]);
	Pdist = P - repmat([1:size(P,1)]',[1 size(P,2)]);
	Pdist(P == 0) = 100;
	Pvis(Pdist == -2) = 1;
	Pvis(Pdist == -1) = 0.5;
	Pvis(cat(3,Z,Pdist == 1,Z)) = 0.5;
	Pvis(cat(3,Z,Pdist == 2,Z)) = 1;
	Pvis(repmat(Pdist == 0,[1 1 3])) = 0.5;
	imshow(Pvis);
	imwrite(Pvis, 'dp_prev.png','Mode','lossless');
end

