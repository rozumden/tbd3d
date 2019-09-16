function [I] = evaluate_vcoeff(isize, tpar, ivl)
% render a trajectory into an image
%
% [I, T] = evaluate_vcoeff(isize, tpar, ivl)
%
% isize ... output image size; [height, width]
% tpar ... trajectory parameters x(t) = a + b*t + c*t^2; format: [a,b,c]
%       could be a cell array of parameters -> then generates multiple
%       curves
% ivl ... interval on which evaluate
%
% H ... rendered trajectory in the image of size isize
%
% Example: 
% tpar = [10 40 0; 40 50 -30];
% [H len] = evaluate_vcoeff([100,100],tpar);


% distance of curve influence
% controls the curve thickness
dthr = 1;
parts = 100;

[~, ~, p1_ext] = postprocc(tpar, [], parts, [0 1], false);
st_ind = max(1, round(parts * ivl(1)));
en_ind = round(parts * ivl(2));
p1 = p1_ext(st_ind:en_ind,:);

[~,~,mask] = renderLineMulti(fliplr(p1), isize);
mask = imdilate(mask, ones(21)); % assume deviation max 10 pixels on each side
[Y X] = find(mask);

% find distances
P = cat(3,X,Y);
D = sqrt(sum((P-permute(p1', [3 2 1])).^2,3));
[m,i] = min(D,[],2);
mask = m<dthr;

% pixel intensity as the length of chord (pixels are assumed circular with radius dthr)
% this is maybe slighty better than the linear model above
T = [X(mask), Y(mask), 2*sqrt(dthr^2 - m(mask).^2)];

I = zeros(isize);
I(sub2ind(isize,T(:,2),T(:,1))) = T(:,3);
