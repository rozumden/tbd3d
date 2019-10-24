function R = createRnMatrix(img_sz, angles)
% creates sparse matrix that averages image over several discrete rotations (used for enforcing rotational symmetry)

if(nargin == 1)
	angles = [16 25 78 152]/180*pi; % +35; selected set of angles to force approx rotational symmetry while more or less avoiding n-fold symmetry (as for [90 120] etc)
end

img_sz = img_sz(1:2);
idx_in = cell(length(angles),1);
idx_out = cell(length(angles),1);

% output coordinates (construct once for all angles)
[idx2 idx1] = meshgrid(1:img_sz(2), 1:img_sz(1));
offset = (img_sz+1)/2; % offset between indices (idx1,2) and coordinates of pixels coordinates for even-sized img are int+.5, for odd-sized are integer
idx = [idx1(:) idx2(:)] - offset;

for i=1:length(angles)
	% no antialiasing, NN 'interpolation'
	% backward rotation - calculate input index for each output index (to avoid vaccuum)

	% rotation matrix
	ca = cos(angles(i)); sa = sin(angles(i));
	R = [ca -sa; sa ca]; % rotates by angle; R.' rotates by -angle

	% input indices (rotated)
	res = round(idx*R + offset); % res contains [in_idx1(:) in_idx2(:)] 1-based coords

	% crop to image dimensions
	keep = all(res >= 1 & res <= img_sz,2);

	% recalculate to linear indices
	idx_out{i} = find(keep);
	idx_in{i} = (res(keep,1)-1)*img_sz(1)+res(keep,2);
end

% 'sum' matrix
R = sparse(cell2mat(idx_out), cell2mat(idx_in), 1/length(angles), prod(img_sz), prod(img_sz));

% averaging
temp = sum(R,2);
w = temp > 0;
R(w,:) = R(w,:)./temp(w);
end
