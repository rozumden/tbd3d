function [pts, I] = trajPredictor(tpar, k, expos, isize)
% render a trajectory into an image
%
% [I, T] = trajRender(isize, tpar, N)
%
% isize ... output image size; [height, width]
% tpar ... trajectory parameters x(t) = a + b*t + c*t^2; format: [a,b,c]
%       could be a cell array of parameters -> then generates multiple
%       curves
% N ... number of points generated on the curve for which we find the
% nearest to every pixel; default 5*<curve length>
%
% I ... rendered trajectory in the image of size isize

% distance of curve influence
% controls the curve thickness
dthr = 1;

if nargin < 3
    expos = 1;
end

if iscell(tpar)
    P = length(tpar);
else
    P = 1;
    tpar = {tpar};
end

% linear->quadratic
w = cellfun(@(x)size(x,2),tpar) < 3;
tpar(w) = cellfun(@(x)[x zeros(size(x,2),1)], tpar(w), 'UniformOutput', false);

% get appropriate parts and proportions
pointsPerUnit = 20; % linearize the curve by this many segments
if k == 0
    inter = [0 1];
else
    t = linspace(0,1,pointsPerUnit);
    t = [ones(length(t),1) t(:) t(:).^2];
    lens = zeros(1,P);
    for c=1:P
        p = t*tpar{c}.'; %%%%% '
        lens(c) = sum(sqrt(sum((p(2:end,:)-p(1:end-1,:)).^2,2)));
    end
    if k > 0
        prop = sum(lens) / lens(end);
        extlen = prop / expos;
        gaplen = (1-expos) * extlen;
        stpoint = 1+gaplen+(k-1)*extlen;
        inter = [stpoint (stpoint+prop)];
        tpar = tpar(end);
    else
        prop = sum(lens) /  lens(1);
        extlen = prop / expos;
        gaplen = (1-expos) * extlen;
        stpoint = gaplen+(abs(k)-1)*extlen;
        inter = [-stpoint -(stpoint+prop)];
        tpar = tpar(1);
    end
    P = 1;
end

% rough curve outline to determine length and relevant points in the output
pts = [];
t = linspace(inter(1),inter(2),pointsPerUnit); 
t = [ones(length(t),1) t(:) t(:).^2];
N = 0; len = 0;
for c=1:P
    p = t*tpar{c}.';
    len = len + sum(sqrt(sum((p(2:end,:)-p(1:end-1,:)).^2,2)));
    pts = [pts; p];
    N(c+1) = round(4*len);
end

if nargin < 4 || nargout < 2
    return
end

[~,~,mask] = renderLineMulti(fliplr(pts), isize);
mask = imdilate(mask, ones(21)); % assume deviation max 10 pixels on each side
[Y X] = find(mask);

% curve parameter
t = zeros(sum(N),1);
C = zeros(2,sum(N));
ind = zeros(sum(N),1);
for c= 1:P
    l = linspace(inter(1),inter(2),N(c+1));
    si = sum(N(1:c));
    t(si+(1:N(c+1)))=l;
    ind(si+(1:N(c+1)))=c;
    % generate N points on the curve
    C(:,si+(1:N(c+1))) = tpar{c}(:,1) + tpar{c}(:,2)*l + tpar{c}(:,3)*l.^2;
end
    
% % coordinates of all pixels in the image
% [X,Y] = meshgrid(1:isize(2),1:isize(1));

% find distances
perm = cat(3,X,Y);
D = sqrt(sum((perm-permute(C, [3 2 1])).^2,3));
% D = sqrt((X(:)-C(1,:)).^2 + (Y(:)-C(2,:)).^2);
[m,i] = min(D,[],2);
mask = m<dthr;

% pixel intensity as the linear function of distance
%T = [X(mask), Y(mask), 1-m(mask)/dthr, l(i(mask)).'];

% pixel intensity as the length of chord (pixels are assumed circular with radius dthr)
% this is maybe slighty better than the linear model above
T = [X(mask), Y(mask), 2*sqrt(dthr^2 - m(mask).^2), t(i(mask)), ind(i(mask))];

I = zeros(isize);
I(sub2ind(isize,T(:,2),T(:,1))) = T(:,3);
end
