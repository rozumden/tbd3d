function [] = make_video_curves(frms, cst, curves, video, sq, folder, gm)
add_original = true;
add_single = false;
add_increasing = true;
add_curves = true;

add_original_slow = true;

add_tsr_dummy = false;
add_removal = true;
add_tsr = true;
add_gt = false;

rsz = 1;
name = sq.short; 
low_contrast = sq.low_contrast;
fps = sq.fps;

fps_use = min([25 fps]);

txt_frms = round(1.2*fps_use);

writer = VID.VideoWriterWrapper('./', name, 'fps', fps_use);
tx = 45;
ty = 45;
font_size = 30;
color_text = [0.6350 0.0780 0.1840];

if add_original
	txt = 'Input Sequence';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(video.get_frame(1), rsz);
	imshow(IM.^gm); hold on
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM.^gm); hold on
		text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
		writer.write();
		clf
	end
end

if add_single
	txt = 'TbD Trajectories';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(video.get_frame(1), rsz);
	imshow(IM.^gm); hold on
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM.^gm); hold on
		if ~isempty(frms{k})
			yellow_at = 0.25;
			b = 1/(1-yellow_at);
			a = -b;
			clr_r = max(0,min(1,a*cst(k)+b));

			a = 1/yellow_at;
			clr_g = max(0,min(1,a*cst(k)));
			clr_b = 0;
	       	clr = [clr_r clr_g clr_b];
	       	frms{k}.show(clr);
		end
		text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
		writer.write();
		clf;
	end	
end

if add_increasing
	txt = 'Accumulated TbD Trajectories';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(video.get_frame(1), rsz);
	imshow(IM.^gm);
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM.^gm);
		for k1 = 1:k
			if ~isempty(frms{k1})
				yellow_at = 0.25;
				b = 1/(1-yellow_at);
				a = -b;
				clr_r = max(0,min(1,a*cst(k1)+b));

				a = 1/yellow_at;
				clr_g = max(0,min(1,a*cst(k1)));
				clr_b = 0;
		       	clr = [clr_r clr_g clr_b];
		       	frms{k1}.show(clr);
			end
		end
		text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
		writer.write();
		clf;
	end
end

if add_curves
	txt = 'Accumulated TbD-NC Trajectories';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(video.get_frame(1), rsz);
	imshow(IM.^gm);
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	xx = {};
	yy = {};
	for k = 1:numel(curves)
		len = 1;
		[y,x] = find(curves(k).Fit);

		if curves(k).fit_iv(1) > 0
			for ck = curves(k).fit_iv(1):curves(k).fit_iv(2)
				IM = imresize(video.get_frame(ck), rsz);
					
				imshow(IM.^gm);
				hold on
				plot(x,y,'.g');
				for kk = 1:numel(xx), plot(xx{kk},yy{kk},'.g'); end
				text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
				writer.write();
				clf;
			end
		else
			imshow(IM.^gm);
			hold on
			plot(x,y,'.g');
			for kk = 1:numel(xx), plot(xx{kk},yy{kk},'.g'); end
			text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
			writer.write();
			clf;
		end
		xx = [xx {x}];
		yy = [yy {y}];
	end
end

if add_original_slow
	txt = 'Input Sequence Slow';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(video.get_frame(1), rsz);
	imshow(IM.^gm);
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM.^gm);
		text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
		for kk = 1:3
			writer.write();
		end
		clf
	end
end

if add_tsr
	n = 240 / fps;

	[~,sname,~] = fileparts(sq.name);
	tdir = load(fullfile(folder, 'templates', [sname '_template.mat']));
	if tdir.scale ~= sq.resize
		tdir.template = imresize(tdir.template, sq.resize / tdir.scale);
	end
	if ~isempty(sq.pad_pcent)
		[template m0] = setupTemplate(tdir.template,sq.pad_pcent);
	else
		[template m0] = setupTemplate(tdir.template);
	end
	r0 = double(size(template,1)/2);
	F = im2double(template);
	M = double(m0);

	out1 = video.get_mat();
	out = VID.remove_fmo(out1,curves,r0); %% remove fmos
	out_removal = VID.slower_video(out,n); 			   %% interpolate without fmos

	if isempty(sq.non_uniform)
		out_tsr = VID.put_fmo(out_removal,curves, M, F, n);  %% interpolate with fmos
	else
		out_tsr = VID.put_fmo_rotation(out_removal,out1,curves, M, F, n);  %% interpolate with fmos and update
	end

	out_tsr(out_tsr > 1) = 1;
	out_tsr(out_tsr < 0) = 0;

	if add_removal
		txt = 'FMO Removal';
		fprintf('%s: %s\n', sq.short, txt);
		IM = imresize(out(:,:,:,1), rsz);
		imshow(IM.^gm);
		text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
		for k = 1:txt_frms, writer.write(); end
		clf;
		for k = 1:size(out, 4)
			IM = imresize(out(:,:,:,k), rsz);
			imshow(IM.^gm);
			text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
			writer.write();
			clf;
		end
	end

	txt = 'Temporal Super-Resolution';
	fprintf('%s: %s\n', sq.short, txt);
	IM = imresize(out_tsr(:,:,:,1), rsz);
	imshow(IM.^gm);
	text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
	for k = 1:txt_frms, writer.write(); end
	clf;
	for k = 1:size(out_tsr, 4)
		IM = imresize(out_tsr(:,:,:,k), rsz);
		imshow(IM.^gm);
		text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
		writer.write();
		clf;
	end
end

if add_gt
	try
		gtvideo = VID.VideoImg([folder '_hs'], sq.short, '');
		n = 240 / fps;
		
		offset = 0;
		if ~isempty(sq.gt_offset), offset = sq.gt_offset; end

		st = (offset + sq.start_ind - 1)*n + 1;
		if isempty(sq.end_ind) 
			en = gtvideo.sz; 
		else
			en = (offset + sq.end_ind - 1) * n;
		end

		txt = 'Ground Truth';
		fprintf('%s: %s\n', sq.short, txt);
		IM = imresize(gtvideo.get_frame(st), sq.resize*rsz);
		imshow(IM.^gm);
		text(tx,4*ty,txt,'FontSize',1.2*font_size,'Color',color_text);
		for k = 1:txt_frms, writer.write(); end
		clf;
		
		
		for k = st:en
			IM = imresize(gtvideo.get_frame(k), sq.resize*rsz);
			imshow(IM.^gm);
			text(tx,ty,txt,'FontSize',font_size,'Color',color_text);
			writer.write();
			clf;
		end
	catch
	end
end

writer.close();

