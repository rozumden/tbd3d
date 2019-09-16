seq(4).low_contrast = 0;
seq(1).no_gt = 1;

figure;

% if isempty(gcp('nocreate')), parpool(numel(inds)); end
% parfor i = 1:numel(seq)
for i = 1:numel(seq)
	warning('off','all');
	disp(seq(i).short);
	gm = 1;
	if seq(i).low_contrast, gm = 1/1.6; end
	t = load(fullfile(folder, seq(i).name));
	if ~isempty(seq(i).end_ind)
		video = imresize(im2double(t.Vk(:,:,:,seq(i).start_ind:seq(i).end_ind)), seq(i).resize);
	else
		video = imresize(im2double(t.Vk(:,:,:,seq(i).start_ind:end)), seq(i).resize);
	end	
	videomat = VID.VideoMat(video);
	seq(i).non_uniform = 1;
	% VIS.make_video(frame{i}, tiou{i}, seq(i).short, video);
	VIS.make_video_curves(frame{i}, tiou{i}, curves{i}, videomat, seq(i), folder, gm);

	fprintf('Done: %s\n', seq(i).short);
end
