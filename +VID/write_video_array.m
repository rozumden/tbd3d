function [] = write_video(video, file, arr)
fps = 5;
format_file = 'Motion JPEG AVI';

v = VideoWriter(file, format_file);
v.FrameRate = fps;

open(v);
for i = 1:size(video,4)
	text_str = ['Value ' num2str(arr(i))];
	RGB = insertText(video(:,:,:,i),[100 100],text_str,'FontSize',18);
	writeVideo(v,RGB);
end
close(v);