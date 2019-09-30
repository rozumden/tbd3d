function img = loadFloFile(fileName)
%LOADFLOFILE Load optical flow in Middlebury benchmark format.

% Author: Damien Teney

fid = fopen(fileName, 'r');
if fid < 0
  img = []; return;
end

tag = fread(fid, 1, 'float32');
width = fread(fid, 1, 'int32');
height = fread(fid, 1, 'int32');

% Sanity checks
tagCheck = 202021.25;
if (tag ~= tagCheck) || (width < 1) || (width > 99999) || (height < 1) || (height > 99999)
  img = []; return;
end

img = single(fread(fid, inf, 'float32'));
fclose(fid);
img = reshape(img, [2, width, height]);
img = permute(img, [3 2 1]) ;

% Remove invalid values
img(img > 999) = NaN ;
