[seq, folder] = EVAL.get_seq(0,'TbD');
tiou_curves = repmat({[]}, 1, numel(seq));
curves_all = repmat({[]}, 1, numel(seq));
rds = repmat({[]}, 1, numel(seq));
if isempty(gcp('nocreate')), parpool(numel(inds)); end
parfor i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).name)
	t = load(fullfile(folder, seq(i).name));
	% EXPER.gen_figure_images(frame{i}, tiou{i}, uint8(imresize(t.Vk(:,:,:,1),scl_fctr)), seq(i).short, scl_fctr/seq(i).resize);
	% VIS.tex_report_image(frame{i}, tiou{i}, uint8(imresize(t.Vk(:,:,:,1),scl_fctr)), seq(i).short, scl_fctr/seq(i).resize);
	im = imresize(t.Vk(:,:,:,seq(i).start_ind+2),seq(i).resize);
	[curves, frms] = FIT.sequence_fit(frame{i}, im);
	if false
		for k = 1:numel(curves)
			% if numel(curves(k).coeff{1}) == 6
			% 	frms_all = [frms{:}];
			% 	rcm = ALG.pixel2cm(mean([frms_all.Radius]),curves(k).coeff{1}(2,3));
			% 	rds{i}(k) = rcm;
			% end
			imshow(im);
			if ~strcmp(curves(k).type,'connect')
				for frm = [frms{ ceil([curves(k).iv(1):curves(k).iv(2)]/8) }]
					frm.show
					% [H,T0,len] = myTrajRender(size2(im), curves(k).coeff, [frm.instance frm.instance+0.1]);
					% [y,x] = find(H);
					% plot(x,y,'.r')
					% plot(frm.Start(1),frm.Start(2),'.b');
					% plot(frm.End(1),frm.End(2),'.b');
				end
			end
			hold on
			plot(curves(k).curve(1,:), curves(k).curve(2,:),'.y');
			[y,x] = find(curves(k).Fit);
			plot(x,y,'.m');
			drawnow
			input('')
		end
	end
	curves_all{i} = curves;
	[tiou_curves{i}] = ALG.gt_cost_iou_curves(curves, frame{i}, par{i}, r0_all(i)); 
	tiou_curves_mean(i) = mean(tiou_curves{i}(3:end));
	disp(tiou_curves_mean(i));
	VIS.tex_report_image_curves(curves_all{i}, uint8(imresize(t.Vk(:,:,:,1),1)), seq(i).short, 1/seq(i).resize);
	% VIS.make_video_curves(frame{i}(1:seq(i).end_frame-1), tiou0(1:seq(i).end_frame-1), curves, seq(i).name, video, seq(i).resize);
end
disp(mean(tiou_curves_mean));
% save('../data/eTbD-T1-LT.mat','tiou_curves','curves_all','rds','tiou_curves_mean','seq','frame','seq','params','cfg')
% save('../data/TbD-T1(0.713)-LT.mat','tiou_curves','curves_all','rds','tiou_curves_mean','seq','frame','seq','params','cfg')

rds = repmat({[]}, 1, numel(seq));
for i = 1:numel(seq)
	for k = 1:numel(curves_all{i})
		if numel(curves_all{i}(k).coeff{1}) >= 6
			frms_all = [frame{i}{:}];
			r = mean([frms_all.Radius]);
			rcm = ALG.pixel2cm(r, norm(curves_all{i}(k).coeff{1}(:,3)), seq(i).fps);
			rds{i}(k) = rcm;
		end
	end
end
