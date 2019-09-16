function [I, T] = trajRender(isize, tpar)
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

if iscell(tpar)
    P = length(tpar);
else
    P = 1;
    tpar = {tpar};
end

% linear->quadratic
w = cellfun(@(x)size(x,2),tpar) < 3;
tpar(w) = cellfun(@(x)[x zeros(size(x,2),1)], tpar(w), 'UniformOutput', false);

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
t = linspace(0,1,20); % linearize the curve by this many segments
t = [ones(length(t),1) t(:) t(:).^2];
N = 0; len = 0;
for c=1:P
    p = t*tpar{c}.';
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
    l = linspace(0,1,N(c+1));
    si = sum(N(1:c));
    t(si+(1:N(c+1)))=l;
    ind(si+(1:N(c+1)))=c;
    % generate N points on the curve
    C(:,si+(1:N(c+1))) = tpar{c}(:,1) + tpar{c}(:,2)*l + tpar{c}(:,3)*l.^2;
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
