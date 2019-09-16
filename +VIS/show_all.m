function [] = show_all(t, seq, frms)
if ~isempty(t)
	image(uint8(imresize(mean(t.Vk,4), seq.resize)));
end
num = 0;
for ki = 1:numel(frms)
	fr = frms(ki);
	for f0 = fr{1}
		% if ~isempty(f0.Speed) && f0.Speed < 0.5 && ki > 1 && ~isempty(frms{ki-1})
			% ALG.update_frame(frame{i}{ki}, frame{i}{ki-1}, 0.9);
		% end
		f0.show;
		f0.instance = ki;
		f0.index = 0;
		% if f0.Preference > 1 && f0.SpeedFull >= 0
		num = num + 1;
		f0.show_number(num);
		f0.index = num;
		% end
	end
end
