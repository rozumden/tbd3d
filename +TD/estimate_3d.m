function [sz] = estimate_3d(Fs,Ms,meth)
if ~exist('meth','var')
	meth = 2;
end
r = size(Ms,1);
sz = zeros(1,size(Fs,4));
rs = [(r/3):0.1:(sqrt(2)*r/2)];
ang = linspace(0, 2*pi, 30);
sang = sin(ang);
cang = cos(ang);

for k = 1:size(Fs,4)
	if size(Ms,4) > 1
		M = Ms(:,:,:,k);
	else
		M = Ms(:,:,k);
	end
	
	if meth == 1
		c = size(M) / 2;
		Vals = [];
		for rd = rs
			Xq = c(1) + 0.5 + rd*sang;
			Yq = c(2) + 0.5 + rd*cang;
			Vq = interp2(M,Xq,Yq,'cubic');
			Vals = [Vals nanmedian( Vq ) ];
		end
		s = find(Vals < 0.5);
		if isempty(s)
			[~,s] = min(Vals);
		end
		s = s(1);
		sz(k) = rs(s);
	elseif meth == 2
		sz(k) = sqrt(sum(M(:))/pi);
	end
end

