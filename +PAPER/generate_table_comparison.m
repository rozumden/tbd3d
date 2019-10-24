function [srctable] = generate_table_comparison()
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

tious = [mean_tiou; mean_tiou_nc; mean_tiou3d; mean_tiou3d_nc; mean_tiou_nc3d3d; mean_tiou3d_nc_oracle; mean_tiou3d_nc3d_oracle; rerr; rerr_est; rerr_gt;];
names = {'TbD', 'TbD-NC', 'TbD', 'TbD-NC', 'TbD-3D', 'TbD-Oracle', 'TbD-3D-Oracle', 'a','a','a'};
which_tiou = logical(ones(numel(names),1));
which_tiou(end-2:end) = 0;

for i = 1:numel(seq)
	[~,seqname,~] = fileparts(seq(i).name);
	C = strsplit(seqname,'_');
	short = C{end};

	c3 = int2str(numel(frame{i}));
	srctable = [srctable short sep c3];

	[max_v] = max(tious(which_tiou,i)); ind = find(max_v == tious(:,i));
	for k = 1:numel(names)
		col = sprintf(['%.' int2str(precision) 'f'],tious(k,i));

		lf = ''; rg = '';
		if any(k == ind)
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
	if mean(tious(k,:)) == max_v_t
		lf = '\\textbf{';
		rg = '}';
	end

	if which_tiou(k)
		meaniou = meaniou(2:end);
	end

	srctable = [srctable sep lf meaniou rg];
end

srctable = [srctable '\\\\\n'];

