function [cs val] = findSignificantRun(s, w, labels, lbl_sum, gap_thresh_small, gap_thresh_big, sig_thresh1, sig_count1, sig_pcent1, sig_thresh2, sig_count2, sig_pcent2)
% note: ma tvrdsi podminky na souvislost small-chunks nez _mini verze, takze se chova trochu jinak, ale s urcitym predzpracovanim funguje podobne (oboje ma vyhody i nevyhody)
% gap_thresh_big - when set to 0 it will only trim the boundary points and remove non-significant chunks on the sides (it is as if gap_thresh_big=Inf but that is not possible, so 0 is used as a flag)
% sig_thresh1 - when set to 0, no significance of chunks is analyzed and simply most valuable run (connected wrt gap_thresh_small) si returned; gap_thresh_big and other params related to significance then dont matter
% 
% cs - returned ordered, so start end end can be determined as cs(1) and cs(end)

% partition into chunks
[s, p] = sort(s); w = w(p); labels = [0; labels(p); 0];
ds = [Inf; s(2:end)-s(1:end-1); Inf]; % ds(i) is spacing before i-th pt and ds(i+1) is spacing after i-th pt
gaps_small = find(ds > gap_thresh_small | labels(2:end)~=labels(1:end-1)); % i-th small-chunk is gaps(i):gaps(i+1)-1
cw = [0; cumsum(w)];
gaps1end = gaps_small(1:end-1); gaps2end = gaps_small(2:end); % auxiliaries

% determine significance of chunks
if(sig_thresh1)
	% chunk is significant is condition(sig_thresh1, sig_count1, sig_pcent1) or  condition(sig_thresh2, sig_count2, sig_pcent2) hold;
	% condition(...1) is "there is at least sig_count1 pixels > sig_thresh1 and all chunk pts form at least sig_pcent1 % of the whole component mass" etc

	% calculate some metrics of mini-chunks
	temp = [0; cumsum(w>=sig_thresh1)]; num1 = temp(gaps2end)-temp(gaps1end); % count of type1 significant pixels for each small-chunk
	temp = [0; cumsum(w>=sig_thresh2)]; num2 = temp(gaps2end)-temp(gaps1end); % count of type2 significant pixels for each small-chunk
	
	if(sig_pcent1 || sig_pcent2)
		temp = labels(gaps1end+1); temp2 = accumarray(temp, cw(gaps2end)-cw(gaps1end), size(lbl_sum)); % sums of pixels of particular labels
		pcent = temp2./lbl_sum; pcent = pcent(temp); % percentage of the sum of the whole component for each small-chunk (lookup via its label)
	else
		pcent = 1;
	end
		
	% state small-chunk significance
	is_sig = (num1 >= sig_count1 & pcent >= sig_pcent1) | (num2 >= sig_count2 & pcent >= sig_pcent2); % is_sig(i) says if i-th small-chunk is significant
	
	% neighboring significant chunks
	temp = [0; find(is_sig); 0]; temp2 = cumsum(is_sig);
	prev_sig = temp(temp2+1); % prev_sig(i) returns idx of closest previous (as in <=) significant segment to the i-th (ie can return the i-th itself); 0 means there is no previous significant chunk
	next_sig = temp(temp2-is_sig+2); % next_sig(i) returns idx of closest next (as in >=) significant segment to the i-th (ie can return the i-th itself); 0 means there is no next significant chunk

	% find big chunks
	if(gap_thresh_big)
		gaps_big = find(ds(gaps_small) > gap_thresh_big); % i-th big chunk runs is all chunks from gaps_big(i) to gaps_big(i+1)-1 (ie indexing is into small chunks, not into pts)
	else % treat the whole set as 1 big chunk (just trim non-sig segments on the sides)
		gaps_big = [1 size(gaps_small,1)]; % all small-chunks
	end
	runs = [next_sig(gaps_big(1:end-1)) prev_sig(gaps_big(2:end)-1)]; % [from to] indices of chunks of significant runs (not pixels)
	runs = runs(all(runs > 0,2), :); % remove invalid runs (not containing any significant chunk)
	runs = [gaps_small(runs(:,1)) gaps_small(runs(:,2)+1)-1]; % convert to px indices ([from to])
else % no significance analysis, simply best run
	runs = [gaps1end gaps2end-1];
end

% calc value of all significant chunks
runs = runs(runs(:,2)-runs(:,1) > 0,:); % remove single-px runs
run_val = cw(runs(:,2)+1)-cw(runs(:,1)); % values of individual significant runs

% choose best
[val, ix] = max(run_val);
cs = p(runs(ix,1):runs(ix,2));
end
