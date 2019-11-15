[seq, folder] = EVAL.get_seq(0,'TbD-3D');
n = 8;
resz = 1;
use_n = false;

% if ~exist('matF','var') || isempty(matF) || any(cellfun(@(x) isempty(x), matF))
	if use_n
		resz = 0.5;
		load('~/projects/data/TbD-3D-n8.mat');
		load('~/projects/data/TbD-3D-n8_matF.mat');
		load('~/projects/data/TbD-3D-n8_post.mat');
	else
		load('~/projects/data/TbD-3D.mat');
	end
% end

tiou = repmat({[]}, 1, numel(seq));
tiou_nc = repmat({[]}, 1, numel(seq));
tiou3d = repmat({[]}, 1, numel(seq));
tiou3d_nc = repmat({[]}, 1, numel(seq));
tiou3d_nc_oracle = repmat({[]}, 1, numel(seq));
tiou3d_nc3d_oracle = repmat({[]}, 1, numel(seq));
tiou_nc3d3d = repmat({[]}, 1, numel(seq));

matF_hs = repmat({[]}, 1, numel(seq));
matM_hs = repmat({[]}, 1, numel(seq));

for i = 1:numel(seq)
	disp(seq(i).name);

	if ~exist('Vk_WB','var') || numel(Vk_WB) < i || isempty(Vk_WB{i}) 
		t = load(fullfile(folder, seq(i).name));
		[Vk{i},Vk_WB{i},PAR{i},V_WB{i}] = generate_lowFPSvideo(t.V,t.POS,t.R,n,resz);
		Vk{i}(:,:,1,:) = 1.9.*Vk{i}(:,:,1,:);
		Vk{i}(:,:,3,:) = 1.8.*Vk{i}(:,:,3,:);
	end

	frms = [frame{i}{:}];
	r_est = max([frms.Radius]);

	[tiou{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i});
	fprintf('[TbD] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou{i}(3:end)]));

	[tiou_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i});
	fprintf('[TbD-NC] Mean TIoU for %s is %.4f\n', seq(i).name, mean([tiou_nc{i}(3:end)]));

	%% TIoU-3D
	[tiou3d{i}] = FIT.gt_cost_iou_3d(frame{i}, PAR{i}, 2);
	fprintf('[TbD] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d{i}(3:end)]));

	[tiou3d_nc{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2);
	fprintf('[TbD-NC] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc{i}(3:end)]));

	%% TIoU-3D TbD-Oracle
	[tiou3d_nc_oracle{i},~,gt_coeffs{i}] = FIT.gt_cost_iou_3d_oracle(r_est, PAR{i});
	fprintf('[TbD-NC-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc_oracle{i}(3:end)]));

	szs_gt{i} = FIT.depth_fit(ind_gt{i}, szs_gt{i});
	[tiou3d_nc3d_oracle{i}] = FIT.gt_cost_iou_3d_oracle(szs_gt{i}, PAR{i},ind_gt{i});
	fprintf('[TbD-NC-3D-Oracle] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou3d_nc3d_oracle{i}(3:end)]));

	%% TbD-3D
	if isempty(curves{i})
		tiou_nc3d3d{i} = zeros(1,numel(frame{i}));
	else
		szs_est{i} = FIT.depth_fit(ind_est{i}, szs_est{i});
		[tiou_nc3d3d{i}] = FIT.gt_cost_iou_3d_curves(curves{i}, frame{i}, PAR{i}, 2, szs_est{i}, ind_est{i});
	end
	fprintf('[TbD-NC-3D] Mean TIoU-3D for %s is %.4f\n', seq(i).name, mean([tiou_nc3d3d{i}(3:end)]));

	[matF_hs{i}, matM_hs{i}] = TD.get_views_hs_3d(im2double(V_WB{i}),curves{i},PAR{i},n,true);
end
EVAL.get_stats;

if false %% reconstruction figure
	for i = [1 4 5]
		PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i});
		PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i}, szs_gt{i}, ind_gt{i});
		PAPER.gen_tex_image_gt(curves{i}, Vk{i}(:,:,:,10), seq(i).name, PAR{i}, szs_est{i}, ind_est{i});
		% PAPER.gen_tex_image_gt(gt_coeffs{i}, Vk{i}(:,:,:,10), seq(i).name, szs_est{i}, ind_est{i});
		% PAPER.generate_imgs(matF_gt{i}, frame{i}, matF_hs{i}, ind_gt{i}(1:8));
	end
	% i=1;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([1:8]));
	% i=7;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([1:8]));

	for i = [1 7]
		scl = 4;
		new_coeff = [gt_coeffs{i}{1}{1}(:,1) gt_coeffs{i}{1}{end}(:,1)];
		new_coeff = {[new_coeff(:,1) new_coeff(:,2)-new_coeff(:,1)]*scl};
		inp = imresize(Vk_WB{i}(:,:,:,1),scl);
		Hall = myTrajRender(size2(inp), new_coeff, [0 1]);
		inp(logical(cat(3,zeros(size(Hall)),Hall>0))) = 255;
		HM = conv2(Hall,ones(scl*size(matM_hs{i}(:,:,:,1))),'same');
		[yy,xx]=find(HM>0);
		for tk = 1:n
			position = new_coeff{1}(:,1) + ((tk)/n)*new_coeff{1}(:,2);
			inp = insertText(inp,position',tk,'BoxOpacity',0,'FontSize',22,'TextColor','white');
			clr = [255 255 255];
			position = round(new_coeff{1}(:,1) + ((tk-0.5)/n)*new_coeff{1}(:,2));
			for ii = -1:1
				for jj = -1:1
					for ui = 1:3, inp(position(2)+ii,position(1)+jj,ui) = clr(ui); end
				end
			end
		end
		inp = inp(min(yy(:)):max(yy(:)),min(xx(:)):max(xx(:)),:);
		imwrite(inp,['~/tmp/recon/in' int2str(i) '.png']);
	end
end

if false %% teaser figure
	tbest = load('../data/TbD-3D_best.mat');
	i=5;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, tbest.frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([31 63 75 171]));
	i=2;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, tbest.frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([1]),{'1use'});
	i=5;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, tbest.frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([31 63]),{'2use','3use'});
	i=4;PAPER.generate_imgs(matF_gt{i}, matM_gt{i}, tbest.frame{i}, matF_hs{i},  matM_hs{i}, ind_gt{i}, ind_gt{i}([3]),{'4use'});
end

if false
	PAPER.print_table_overleaf(PAPER.generate_table_comparison());
	PAPER.print_table_overleaf(PAPER.generate_table_comparison_v2());
end

if false %% new teaser
	resz = 0.5;
	n = 2;
	load('~/projects/data/TbD-3D-n2.mat');
	load('~/projects/data/TbD-3D-n2_matF.mat');
	load('~/projects/data/TbD-3D-n2_post.mat');
	i = 5;
	% imshow(imresize(Vk_WB{i}(:,:,:,1),resz));
	% VIS.show_all([], seq(i), frame{i});
	figure; 
	image(ones(size(imresize(Vk_WB{i}(:,:,:,1),resz))));
	hold on;
	VIS.curve(curves{i});
	pbaspect([size(Vk{i},2)/size(Vk{i},1) 1 1]);
	set(gca,'linewidth',3);
	xticks([0:50:size(Vk{i},1)]);
	yticks([0:50:size(Vk{i},2)]);
	set(gca,'FontSize',15);
	saveas(gcf,'~/tmp/teaser_traj.png');
	% VIS.curve(gt_coeffs{i});

	% VIS.curve3d(curves{i},szs_gt{i},ind_gt{i});
	VIS.curve6d(curves{i},szs_gt{i},ind_gt{i},matF_gt{i});
	pbaspect([size(Vk{i},2)/size(Vk{i},1) 1 10]);
	set(gca,'linewidth',3);
	xticks([0:50:size(Vk{i},1)]);
	yticks([0:50:size(Vk{i},2)]);
	set(gca,'FontSize',15);
	box on
	ax = gca;
	ax.BoxStyle = 'full';
	axis tight;
	axis equal
    view([-50,-45]);
	saveas(gcf,'~/tmp/teaser_3dtraj.png');


end