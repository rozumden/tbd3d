function [srctable] = generate_table_comparison_v2()
load('~/projects/data/TbD-3D.mat');

precision = 3;

sep = ' & ';
srctable = '';

[seq, folder] = EVAL.get_seq(0,'TbD-3D');
n = 8;
resz = 1;
Vk = repmat({[]}, 1, numel(seq));
Vk_WB = repmat({[]}, 1, numel(seq));
PAR = repmat({[]}, 1, numel(seq));
V_WB = repmat({[]}, 1, numel(seq));
Ms = repmat({[]}, 1, numel(seq));
EVAL.run_tbd3d_soft;

disp('--------------------------------------------------');

EVAL.get_stats;

if false
	% load('../data/TbD-3D-n1.mat');
	load('../data/TbD-3D-n1_matF.mat');
	load('../data/TbD-3D-n1_post.mat');
	EVAL.run_tbd3d_rotation;
	EVAL.run_fps_rotation;
else
	load('../data/TbD-3D-rotations-abs.mat');
end
err_or_u = rot_averages.alla(1).err_or_u;
err_u = rot_averages.alla(1).err_u;
err_or_rot = rad2deg(rot_averages.alla(1).err_or_rot);
err_rot = rad2deg(rot_averages.alla(1).err_rot);

do_bold = true;

% mean_tiou3d_nc_oracle
tious = [mean_tiou3d; mean_tiou3d_nc; mean_tiou_nc3d3d; ...
		mean_tiou3d_nc3d_oracle; rerr; rerr_est; rerr_gt; ...
		err_u; err_or_u; err_rot; err_or_rot];
names = {'TbD', 'TbD-NC', 'TbD-3D', 'TbD-3D-Oracle', 'a','a','a','a','a','a','a'};
which_tiou = logical(ones(numel(names),1));
which_tiou(5:end) = 0;

bpnts = [1 5 8 10 12];
where_max = [1 0 0 0];

for i = 1:numel(seq)
	[~,seqname,~] = fileparts(seq(i).name);
	C = strsplit(seqname,'_');
	short = C{end};

	c3 = int2str(numel(frame{i}));
	srctable = [srctable short sep c3];

	ind = [];
	ki = 1;
	for k = 1:numel(names)
		if any(k == bpnts)
			if where_max(ki)
				[max_v] = max(tious(bpnts(ki):(bpnts(ki+1)-1),i)); 
			else
				[max_v] = min(tious(bpnts(ki):(bpnts(ki+1)-1),i)); 
			end
			ind = find(max_v == tious(:,i));
			ki = ki + 1;
		end
		col = sprintf(['%.' int2str(precision) 'f'],tious(k,i));

		lf = ''; rg = '';
		if do_bold &&  any(k == ind)
			lf = '\\textbf{';
			rg = '}';
		end

		if which_tiou(k)
			col = col(2:end);
		end
		srctable = [srctable sep lf col rg];
	end
	srctable = [srctable '\\\\  \n '];
end

meanfrm = sprintf('%.0f',mean(arrayfun(@(x) numel(x{1}), frame)));
srctable = [srctable '\\hline \n  Average ' sep meanfrm];

mean_ttt = mean(tious');
[max_v_t] = max(mean_ttt(which_tiou));
for k = 1:numel(names)
	meaniou = sprintf(['%.' int2str(precision) 'f'],mean(tious(k,:)));
	lf = ''; rg = '';
	if false && mean(tious(k,:)) == max_v_t
		lf = '\\textbf{';
		rg = '}';
	end

	if which_tiou(k)
		meaniou = meaniou(2:end);
	end

	srctable = [srctable sep lf meaniou rg];
end

srctable = [srctable '\\\\\n'];

