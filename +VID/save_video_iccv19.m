inds = [1:numel(seq)];
frame{5}{10} = []; frame{5}{11} = []; 
frame{4}{36} = []; frame{4}{37} = []; frame{4}{38} = []; frame{4}{31} = [];
frame{9}{8} = []; 
gm = 1/1.6
for i = inds
	if i > 3, gm = 1; end
	t = load(fullfile(folder, seq(i).name));
	if ~isempty(seq(i).end_ind)
		video = VID.VideoMat(imresize(im2double(t.Vk(:,:,:,seq(i).start_ind:seq(i).end_ind)).^gm, seq(i).resize));
	else
		video = VID.VideoMat(imresize(im2double(t.Vk(:,:,:,seq(i).start_ind:end)).^gm, seq(i).resize));
	end	
	VIS.make_video(frame{i}, tiou{i}, seq(i).short, video);
end