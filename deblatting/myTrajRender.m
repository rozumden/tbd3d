function [I, T, len] = trajRender(isize, tpar, ivl)
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
% T ... matrix of pixels on the trajectory with weights; format: [x,y,w,t,ind] 
%
% Example: 
% tpar = [10 40 0; 40 50 -30];
% [I T] = trajRender([100,100],tpar);


% distance of curve influence
% controls the curve thickness
dthr = 1;

if nargin < 3
    ivl = [0 1];
end
ivl_len = ivl(2) - ivl(1);

if iscell(tpar)
    P = length(tpar);
else
    P = 1;
    tpar = {tpar};
end

% linear->quadratic
% dof = 6;
% w = cellfun(@(x)size(x,2),tpar) < dof;
% tpar(w) = cellfun(@(x)[x zeros(2,dof-size(x,2))], tpar(w), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % determine the number of points on curve, if N is not defined
% if nargin<3 ||  isempty(N)
%     dt = 1/100;
%     N = zeros(1,P+1);
%     for c = 1:P
%         % curve length in pixels
%         d = sum(sqrt(sum((tpar{c}(:,2)+tpar{c}(:,3)*[0:dt:1]).^2))*dt);
%         % generate 5x more points than is the length; probably too much!!
%         if (d<1)
%             d = 1;
%         end
%         N(c+1) = round(5*d);
%     end
% else
%     if size(N)==1
%         N = [0, repmat(N,1,P)];
%     end
% end

% rough curve outline to determine length and relevant points in the output
pts = [];
N = 0; len = 0;
n_seg = 20;
if ivl_len < 1
    n_seg = 100;
end

for c=1:P
    max_power = (numel(tpar{c}) / 2) - 1;
    ti = linspace(ivl(1),ivl(2),n_seg*ivl_len); % linearize the curve by this many segments
    t = [ones(length(ti),1)];
    for pow1 = 1:max_power
        t = [t ti(:).^pow1];
    end

    p = t*tpar{c}.'; %%%% '
    len = len + sum(sqrt(sum((p(2:end,:)-p(1:end-1,:)).^2,2)));
    pts = [pts; p];
    N(c+1) = round(4*len);
end
[~,~,mask] = renderLineMulti(fliplr(pts), isize);
mask = imdilate(mask, ones(21)); % assume deviation max 10 pixels on each side
[Y X] = find(mask);

% curve parameter
t = zeros(sum(N),1);
C = zeros(2,sum(N));
ind = zeros(sum(N),1);
for c= 1:P
    l = linspace(ivl(1),ivl(2),N(c+1));
    si = sum(N(1:c));
    t(si+(1:N(c+1)))=l;
    ind(si+(1:N(c+1)))=c;
    % generate N points on the curve
    max_power = (numel(tpar{c}) / 2) - 1;
    C(:,si+(1:N(c+1))) = tpar{c}(:,1)*ones(size(l));
    for pow1 = 1:max_power
        C(:,si+(1:N(c+1))) = C(:,si+(1:N(c+1))) + tpar{c}(:,pow1+1)*l.^pow1;
    end
end
    

% % coordinates of all pixels in the image
% [X,Y] = meshgrid(1:isize(2),1:isize(1));

% find distances
P = cat(3,X,Y);
D = sqrt(sum((P-permute(C, [3 2 1])).^2,3));
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
