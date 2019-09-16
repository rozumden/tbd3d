function [] = make_video(frms, cst, name, video, rsz)
if nargin < 5
	rsz = 1;
end
add_original = true;
add_single = true;
add_increasing = true;

writer = VID.VideoWriterWrapper('~/projects/vis/', name, 'fps', 10);

if add_original
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM);
		writer.write();
		clf
	end
end

if add_single
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM);
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

		writer.write();
		clf;
	end	
end

if add_increasing
	clf;
	for k = 1:numel(frms)
		IM = imresize(video.get_frame(k), rsz);
		imshow(IM);
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
		writer.write();
		clf;
	end
end

writer.close();

