
for i = 1:numel(seq)
	clf;
	imshow(Vk{i}(:,:,:,10));
	hold on;
	VIS.show_all([], seq(i), frame{i});
	VIS.curve(curves{i});
	pause
	VIS.curve(gt_coeffs{i});
	pause
end