function [] = write_video(video,file, format_file, fps)
if nargin < 4
	fps = 30;
end

if nargin < 3
	format_file = 'Motion JPEG AVI';
end

v = VideoWriter(file,format_file);
% v.LosslessCompression = true;
v.FrameRate = fps;
v.Quality = 100;
% v.VideoCompressionMethod = 'none';

open(v);
for i = 1:size(video,4)
	writeVideo(v,video(:,:,:,i));
end
close(v);
