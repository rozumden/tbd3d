function [matFcell] = matF_rescale(matF, szs)
matFcell = cell(1,size(matF,4));
radmax = size(matF,1)/2;
for kk = 1:size(matF,4)
	rad = ceil(szs(kk));
	if rad < radmax
		matFcell{kk} = matF(radmax-rad:radmax+rad,radmax-rad:radmax+rad,:,kk);
	else
		matFcell{kk} = matF(:,:,:,kk);
	end
end