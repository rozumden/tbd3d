seq = [];
seq(numel(seq)+1).name = 'ball3d_synth_01.mat';	
seq(numel(seq)+1).name = 'ball3d_synth_02.mat';	
seq(numel(seq)+1).name = 'ball3d_synth_03.mat';	
seq(numel(seq)+1).name = 'ball3d_img_01.mat';	
seq(numel(seq)+1).name = 'ball3d_img_02.mat';	
seq(numel(seq)+1).name = 'ball3d_img_03.mat';	
seq(numel(seq)+1).name = 'ball3d_img_04.mat';	
folder = '/mnt/lascar/rozumden/dataset/synth';

[params, cfg] = EVAL.get_params();
[seq, folder] = EVAL.get_seq(0,'synth');

n = 8;

% for i = 1:numel(seq)
for i = 5
	disp(seq(i).name);
	t = load(fullfile(folder, seq(i).name));
	m = double(diskMask([],t.radius));
	[Fs,Ms,HR,roi] = estimateFM_hierarchy(t.img, t.bg, { fliplr(t.coeffs)'}, n, m, []);
	% h = evaluate_vcoeff(size2(t.img), { fliplr(t.coeffs)'}, [0 1]);
	
	% montage(Fs, 'Size', repmat(ceil(sqrt(size(Fs,4))),[1 2]));
	montage(HR, 'Size', [size(HR,4)/n n]);
	
	matF = zeros([size(m) 3 10]);
	for k = 1:size(matF,4)
		[y,x] = find(sum(abs(t.snapshots{k}),3) > 0);
		rect = [min(x)+1 min(y)+1 2*t.radius-1 2*t.radius-1];
		matF(:,:,:,k) = imcrop(t.snapshots{k}, rect);
	end
	Fs0 = matF(:,:,:,([1:2 4:6 8:10]));
	montage(cat(4,Fs0,Fs), 'Size', [2 size(Fs,4)]);
	%res = TD.reconstruct(matF,true);
end	

