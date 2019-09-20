%% STEP 1: create b, fg, m - rotating object (Example 1)
b = linscale(double(imread('bg.jpg'))); % background
b = b(1:256,1:256,:);
% number of segments
NS = 20;
TF = linscale(double(imread('simpleF.png')));
TM = linscale(double(imread('simpleM.png')));
f = cell(1,NS);
m = cell(1,NS);
a = linspace(0,30,NS);
% generate NS rotated samples of fg
for i = 1:NS
    f{i} = imrotate(TF,a(i),'bilinear','crop');
    m{i} = TM; % same mask
    % m{i} = imrotate(TM,a(i),'bilinear','crop');
end

%% STEP 1: create b, fg, m - object with changing size (Example 2)
b = linscale(double(imread('bg.jpg'))); % background
b = b(1:256,1:256,:);
% number of segments
NS = 20;
TF = linscale(double(imread('simpleF.png')));
TM = linscale(double(imread('simpleM.png')));
f = cell(1,NS);
m = cell(1,NS);
a = linspace(1,0.6,NS);
RA = imref2d([size(TF,1),size(TF,2)],[-1 1],[-1 1]);
% generate NS rotated samples of fg
for i = 1:NS
    T = affine2d([a(i) 0 0; 0 a(i) 0; 0 0 1]);
    f{i} = imwarp(TF,RA,T,'bicubic','OutputView',RA);
    m{i} = imwarp(TM,RA,T,'bicubic','OutputView',RA);
end

%% STEP 2: create blurred image
fdelta = zeros(size(f{1}));
fdelta(floor(size(f{1},1)/2)+1, floor(size(f{1},2)/2)+1,:, :) = 1;
N = 100; % sampling of the blur curve
speed = 10; % length of the blur
% blur parameters are functions of time: angle(t), [x(t), y(t), z(t)],
% depth, position in the image (0,0) is the image center, image size


% generate NS segments of h
par = [linspace(0,3*speed,N); linspace(0,2*speed,N)+10; linspace(0,4*speed,N)+50; linspace(0,0,N)];
h = zeros(size(b,1),size(b,2),NS);
for i = 1:NS
    range = round((i-1)*N/NS)+1:round(i*N/NS);
    %range = round((i-1)*N/NS)+1;
    h(:,:,i) = genstreak(par(1,range), par(2:4,range), 1, [0 0].', [size(b,1), size(b,2)]);
end
h = h/NS;
% blurred foreground
snr = 50;
F = 0;
for i = 1:NS
    F = F + fft2(h(:,:,i),size(b,1),size(b,2)).*fft2(f{i},size(b,1),size(b,2));
end
%vn = sqrt(var(F(:),1)/(10^(snr/10)));
%sigma2 = vn^2;
%F = F + vn*randn(size(F));
% blurred mask
M = 0;
for i = 1:NS
    M = M + fft2(h(:,:,i),size(b,1),size(b,2)).*fft2(m{i},size(b,1),size(b,2));
end

% final image
g = (1-real(ifft2(M))).*b + real(ifft2(F));
vn = sqrt(var(g(:),1)/(10^(snr/10)));
sigma2 = vn^2;
g = g + vn*randn(size(g));
bp = g;

% shift h by size(M)/2, because fft2 was used insted of conv to generate g 
h = circshift(h,floor([size(m{1},1)/2, size(m{1},2)/2]));
%% STEP 3: run deblatting
sM = size(f{1}) + [3 3 0];
%hinit = repmat(ones(size(h,1),size(h,2))/(size(h,1)*size(h,2)),1,1,size(h,3));
R = 4;
NSR = NS/R;
hinit = zeros(size(h,1),size(h,2),NSR);
for i = 1:NSR
    hinit(:,:,i) = sum(h(:,:,[1:R]+(i-1)*R),3);
end
finit = repmat({ones(sM)},1,NSR); % initial segments of H; no update of H in convsparseCG_v1z
%finit = f;

minit = repmat({0},1,NSR); % independent masks
%minit = {0}; % force single mask

[hr,fr,mr,report] = convsparseCG_v1z(g,b,finit,minit,hinit);

%% display results
mrp = repmat(mr,1,1,1,3);
frp = fr; frp(mrp>0.4) = fr(mrp>0.4)./mrp(mrp>0.4);
figure; subplot(211); dispIm(reshape(frp,size(frp,1),size(frp,2)*NSR,[]));
subplot(212); dispIm(reshape(mr,size(mr,1),size(mr,2)*NSR,[]));

% [Fs,Ms,roi] = estimateFM_motion_template_pw(g, b, hinit, [], zeros([size2(finit{1}) 1 NSR]), [], [], ...
        % 'alpha', 3^-10, 'alpha_m', 3^-9, 'lambda', 1e-3, 'maxiter', 20, 'cg_maxiter', 50, 'cg_tol', 1e-6);
