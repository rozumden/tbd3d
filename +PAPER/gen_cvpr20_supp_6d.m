function [] = gen_cvpr20_supp_6d(curves, szs_gt, ind_gt, matF_gt, Vk)
resz = 0.5;
n = 2;
load('~/projects/data/TbD-3D-n2.mat');
load('~/projects/data/TbD-3D-n2_matF.mat');
load('~/projects/data/TbD-3D-n2_post.mat');
i = 5;



VIS.curve6d_video(curves{i},szs_gt{i},ind_gt{i},matF_gt{i},Vk{i},resz);
