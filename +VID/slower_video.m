function out = slower_video(video,n)
out = average(video,n);

% if isempty(dets)
% 	out = average(video,n);
% elseif isempty(n)
% 	out = delete_object(video,dets);
% else
% 	out = cut_fast(video,dets,n);
% end


function out = average(video,n)
previous = video(:,:,:,1);
out = previous;
for i = 2:size(video,4)
	% disp(int2str(i));
	current = video(:,:,:,i);
	for j = 1:(n-1)
		out(:,:,:,(i-2)*n + j) = ((n-j)/n)*previous + (j/n)*current;
	end
	out(:,:,:,(i-1)*n) = current;
	previous = current;
end

function out = delete_object(video,dets)
frame0 = video(:,:,:,1);
out = frame0;
se = strel('square',11);
for i = 2:size(video,4)
	frame = video(:,:,:,i);
	if ~isempty(dets{i}) && ~isempty(dets{i}.d)
		BW = logical(zeros(size(video,1),size(video,2)));
		BW(dets{i}.d.AllPixels) = 1;
		BW2 = imdilate(BW,se);
		mask = repmat(BW2,[1 1 3]);
		frame(find(mask)) = frame0(find(mask)); 
	end
	out(:,:,:,i) = frame;
	frame0 = frame;
end

function out = cut_fast(video,dets,n)
back = video(:,:,:,1);
se = strel('square',11);
if ~isempty(dets{1}) && ~isempty(dets{1}.d)
	BW = logical(zeros(size(video,1),size(video,2)));
	BW(dets{1}.d.PixelIdxList) = 1;
	BW2 = imdilate(BW,se);
	mask = repmat(BW2,[1 1 3]);
	front = video(:,:,:,2);
	back(find(mask)) = front(find(mask));
end
out = back;
last = [];
for i = 2:size(video,4)
	disp(int2str(i));
	front = video(:,:,:,i);
	front_show = front;
	ball = [];
	if ~isempty(dets{i}) && ~isempty(dets{i}.d)
		offset = size(video,1)*size(video,2);
		BW = logical(zeros(size(video,1),size(video,2)));
		BW(dets{i}.d.PixelIdxList) = 1;
		BW2 = imdilate(BW,se);
		mask = repmat(BW2,[1 1 3]);
		im = im2double(back); im(~mask) = 0;
		bg = im2double(front); bg(~mask) = 0;
		dif = im - bg;
		% for a = 0.2
		% 	obj = dif/a + bg;
		% 	color = mean(reshape(obj(mask),[],3));
		% 	color = min([1 1 1],color);
		% 	color = max([0 0 0],color);
			color = [1 1 1];
		% end
		front(find(mask)) = back(find(mask));
		% alpha = min(1,a*(n-1));
		alpha = 0.6;
		if isempty(last)
			last = dets{i}.d.PixelList;
		else
			[x,y] = ind2sub(size(BW),dets{i}.d.Trajectory);
			pixel = [y x]';
			d = pdist2(last',pixel','euclidean','Smallest',1);
			all_pixels = dets{i}.d.PixelIdxList;
			if isfield(dets{i}.d,'inter_traj') && ~isempty(dets{i}.d.inter_traj)
				pixel = [pixel dets{i}.d.inter_traj.pixels'];
				p = dets{i}.d.inter_traj.extended;
				query = sub2ind(size(BW),p(:,2),p(:,1));
				all_pixels = [all_pixels query'];
				if i < size(video,4) && ~isempty(dets{i+1}) && ~isempty(dets{i+1}.d)
					all_pixels = [all_pixels dets{i+1}.d.PixelIdxList];
				end
			end
			[v1,id1] = min(d);
			[v2,id2] = max(d);
			first = pixel(:,id1);
			last = pixel(:,id2);
			if v1 < 100
				piece = floor(size(pixel,1)/n);
				for k = 1:(n-1)
					[~,idx] = pdist2(pixel',first','euclidean','Smallest',piece+1);
					traj = pixel(:,idx(1:(end-1)));
					first = pixel(:,idx(end));
					inl = logical(ones(1,size(pixel,1)));
					inl(idx(1:(end-1)))= 0;
					pixel = pixel(:,inl);
					mask = gen_ball(traj,dets{i}.d.Radius/2,all_pixels,size(BW));
					in = find(mask > 0);
					ball{k}.in = in;
					ball{k}.values = mask(in);
				end
				mask = gen_ball(pixel,dets{i}.d.Radius/2,all_pixels,size(BW));
				in = find(mask > 0);
				ball{n}.in = in;
				ball{n}.values = mask(in);
			end
		end
	end
	
	for j = 1:n
		if ~isempty(ball) && ~isempty(ball{j}.in)
			in = ball{j}.in;
			values = ball{j}.values;
			object = uint8(zeros(size(video,1),size(video,2),3));
			colors = 255*repmat(color,1,numel(in)/3);
			front_show = double(front);
			back_show = back;
			object(in) = uint8((1-alpha*values).*front_show(in) + alpha*values.*colors(:));
			front_show(in) = 0;
			back_show(in) = 0;
			out(:,:,:,(i-2)*n + j) = ((n-j)/n)*back_show + (j/n)*uint8(front_show) + object;
		else
			out(:,:,:,(i-2)*n + j) = ((n-j)/n)*back + (j/n)*front;
		end
	end
	back = front;
end


function mask = gen_ball(traj,r,idx,sz)
mask = [];
if isempty(traj)
	return
end
[y,x] = ind2sub(sz,idx);
query = [x' y'];
mask = double(zeros(sz));
d = pdist2(traj,query,'euclidean','Smallest',1);
dn = d/r;
dn(dn > 1) = 1;
vals = (1-dn).^.5;
mask(idx) = vals;
mask = repmat(mask,[1 1 3]);