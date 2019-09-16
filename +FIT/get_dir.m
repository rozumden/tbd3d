function [direc, speed] = get_dir(cf0, cf1, bb0, bb1, f0, f)
coe0 = cellfun(@(x)fliplr(x), cf0, 'UniformOutput', false);
coe1 = cellfun(@(x)fliplr(x), cf1, 'UniformOutput', false);
if numel(bb0) == 4
	p0_st = coe0{1}(1,:) + [bb0(1:2)] - 1;
	p0_end = sum(coe0{end}) + [bb0(1:2)] - 1;
	
else
	p0_st = TRACK.change_coor(bb0, coe0{1}(1,2), coe0{1}(1,1))';
	temp = sum(coe0{end});
	p0_end = TRACK.change_coor(bb0, temp(2), temp(1))';
end

if numel(bb1) == 4
	p1_st = coe1{1}(1,:) + [bb1(1:2)] - 1;
	p1_end = sum(coe1{end}) + [bb1(1:2)] - 1;
else
	p1_st = TRACK.change_coor(bb1, coe1{1}(1,2), coe1{1}(1,1))';
	temp = sum(coe1{end});
	p1_end = TRACK.change_coor(bb1, temp(2), temp(1))';
end

points = [p0_st; p0_end];

dist_st = sqrt(sum( ((p1_st - points).^2)' ));
dist_end = sqrt(sum( ((p1_end - points).^2)' ));
f0.Direction = 1;
if min(dist_st) < min(dist_end)
	direc = 1;
	[speed, ind] = min(dist_st);
	if nargin == 6
		f0.Start = points(mod(ind,2)+1,:);
		f0.End = points(ind,:);
		if ind == 1
			f0.Direction = -1;
		end
		f.Start = p1_st;
		f.End = p1_end;
	end
else
	direc = -1;
	[speed, ind] = min(dist_end);
	if nargin == 6
		f0.Start = points(mod(ind,2)+1,:);
		f0.End = points(ind,:);
		if ind == 1
			f0.Direction = -1;
		end
		f.Start = p1_end;
		f.End = p1_st;
	end
end

if nargin == 6
	f0.Centroid = (f0.Start + f0.End) / 2;
	f.Centroid = (f.Start + f.End) / 2;
end

