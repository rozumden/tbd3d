% DATASET
%
% Camera setup:
% Basler resolution: 1920x800 (centerY . . offsetY = 208)
% 245fps
% exposure time: 4020us
% RGB8
% no light source (completely raw)
%
% Videos:
% 1) objects thrown (kicked) perpendicularly to the optical axis 
% out1, out2: volleyball
% outa1: aerobie (strong deformations during bounces)
% outb1: black soft ball
% outf1: football 
% 
% 2) objects thrown (kicked) arbitrarily (strong depth effect)
% depthf1, depthf2, depthf3: football 
% (depthf2 – incorrect detection around the contact with the wall)
% depth2: volleyball
% depthb2: black soft ball

%%  depthf1 video football 
% OK
clear dataset
dataset(1).fn = 'depthf1';
dataset(1).startindex = 51;
dataset(1).endindex = 420;
dataset(1).sensitivity = [0.90 0.95];
dataset(1).radii = [18 32]; %[10 16];
dataset(1).sensitivity2 = [0.98 0.99];
dataset(1).smooth = 155;
% depthf2 video football 
% slightly incorrect detection during bounce
dataset(end+1).fn = 'depthf2';
dataset(end).startindex = 1;
dataset(end).endindex = 400;
dataset(end).sensitivity = [0.90 0.95];
dataset(end).radii = [18 32];
dataset(end).sensitivity2 = [0.98 0.99];
dataset(end).smooth = 155;
% depthf3 video football 
% last 68 frames removed (ball size detection failed)
dataset(end+1).fn = 'depthf3';
dataset(end).startindex = 1;
dataset(end).endindex = 368-68;
dataset(end).sensitivity = [0.90 0.95];
dataset(end).radii = [18 42];
dataset(end).sensitivity2 = [0.98 0.99];
dataset(end).smooth = 155;
% depth2 video volleyball 
% Note: first 20 and last 80 frames were removed -- incorrect ball size detected
% OK
dataset(end+1).fn = 'depth2';
dataset(end).startindex = 61+20;
dataset(end).endindex = 550-80;
dataset(end).sensitivity = [0.90 0.95];
dataset(end).radii = [20 40]; %[10 16];
dataset(end).sensitivity2 = [0.98 0.99];
dataset(end).smooth = 105;
% depth2 video black ball
% Some frames at the end are without any detection 
dataset(end+1).fn = 'depthb2';
dataset(end).startindex = 1;
dataset(end).endindex = 650;
dataset(end).sensitivity = [0.96 0.98];
dataset(end).radii = [18 32]; %[10 16];
dataset(end).sensitivity2 = [0.98 0.99];
dataset(end).smooth = 105;
%  out1 video volleyball 
% OK
dataset(end+1).fn = 'out1';
dataset(end).startindex = 18;
dataset(end).endindex = 475;
dataset(end).sensitivity = [0.95 0.97]; %[0.9 0.95];
dataset(end).radii = [25 30]; %[10 15];
dataset(end).sensitivity2 = [0.98 0.99];
% out2 video volleyball 
% OK
dataset(end+1).fn = 'out2';
dataset(end).startindex = 81;
dataset(end).endindex = 480;
dataset(end).sensitivity = [0.95 0.97]; %[0.98 0.99]; %[0.9 0.95];
dataset(end).radii = [25 30]; %28; %[10 15];
dataset(end).sensitivity2 = [0.98 0.99];
%  out3 video volleyball % NOT TESTED
%dataset(end+1).fn = 'out3';
%dataset(end).startindex = 1;
%dataset(end).endindex = 500;
%  video aerobie 
% OK
dataset(end+1).fn = 'outa1'; 
dataset(end).startindex = 1;
dataset(end).endindex = 380;
dataset(end).sensitivity = [0.99 0.995];
dataset(end).radii = 36; %[15 20];
%  video black ball 
% OK
dataset(end+1).fn = 'outb1';
dataset(end).startindex = 15;
dataset(end).endindex = 344;
dataset(end).sensitivity = [0.97 0.98];
dataset(end).radii = [18 23]; %21; %[7 12];
dataset(end).sensitivity2 = [0.98 0.99];
%  video football 
% OK
dataset(end+1).fn = 'outf1';
dataset(end).startindex = 21;
dataset(end).endindex = 500;
dataset(end).sensitivity = [0.94 0.97];
dataset(end).radii = [25 30];
dataset(end).sensitivity2 = [0.98 0.99];

dataset_index = find(ismember({dataset.fn},'depthf3'));
index = dataset_index;
fd = dataset(index);
%%
for index = dataset_index
fd = dataset(index);
disp(fd.fn);
I = readRAWBasler(['~/tmp/3D_bounce_test/',fd.fn,'.dat']);
I = I(:,:,fd.startindex:fd.endindex);
N = size(I,3);
DSfactor = 1;
V = zeros(size(I,1)/2/DSfactor, size(I,2)/2/DSfactor, 3, N,'uint8');
for i = 1:N
    T = double(I(:,:,i));
    T = cat(3,T(1:2:end,1:2:end), ...
        0.5*(T(1:2:end,2:2:end)+T(2:2:end,1:2:end)), T(2:2:end,2:2:end));
    T = convn(T,ones(DSfactor)/DSfactor^2,'valid');
    V(:,:,:,i) = uint8(T(1:DSfactor:end,1:DSfactor:end,:));
end
B = median(V,4);

% white balance & gamma correction
WB = [2 1 2]; gamma_coef = 0.4;


%% 1st stage: find object position and radius using Hough transform
POS = zeros(2,N);
iR = zeros(1,N);
radii = fd.radii;
sensitivity = fd.sensitivity;
for i = 1:N
    % L1 norm & gamma correction
    G = linscale(sum(abs(double(V(:,:,:,i))-double(B)),3)).^gamma_coef;
    % Circ-Hough transform
    [c,r,m] = imfindcircles(G,radii,'Sensitivity', sensitivity(1));
    disp(num2str(i));
    if isempty(m)
        disp('Trying higher sensitivity');
        [c,r,m] = imfindcircles(G,radii,'Sensitivity',sensitivity(2));
    end
    if length(m) ~= 1
        disp(['No. of detections: ',num2str(length(m))]);
    end
    if length(m)>1
        figure(1); dispIm(G); hold on; viscircles(c, r,'EdgeColor','b');
        %keyboard;
    end
    if isempty(m)
        warning('No detection!');
        POS(:,i) = nan;
        iR(i) = nan;
    else
        POS(:,i) = c(1,:);
        iR(i) = r(1);
    end
end
figure(2); clf; plot(POS(1,:),POS(2,:),'x-'); set(gca, 'YDir','reverse');
R = iR;


%% 2nd stage: smoothen radii and re-run Hough detection but this time with fixed radii
if ~isempty(fd.sensitivity2)
if isempty(fd.smooth)
    smoothSize = 155;
else
    smoothSize = fd.smooth;
end
R = smooth(iR,smoothSize,'rlowess');
POS = zeros(2,N);
sensitivity2 = fd.sensitivity2;
for i = 1:N
    % L1 norm & gamma correction
    G = linscale(sum(abs(double(V(:,:,:,i))-double(B)),3)).^gamma_coef;
    % Circ-Hough transform
    [c,r,m] = imfindcircles(G,R(i),'Sensitivity', sensitivity2(1));
    disp(num2str(i));
    if isempty(m)
        disp('Trying higher sensitivity');
        [c,r,m] = imfindcircles(G,R(i),'Sensitivity',sensitivity2(2));
    end
    if length(m) ~= 1
        disp(['No. of detections: ',num2str(length(m))]);
    end
    if length(m)>1
        figure(1); dispIm(G); hold on; viscircles(c, r,'EdgeColor','b');
        %keyboard;
    end
    if isempty(m)
        warning('No detection!');
        POS(:,i) = nan;
    else
        POS(:,i) = c(1,:);
    end
end
figure(2); hold on; plot(POS(1,:),POS(2,:),'rx-'); set(gca, 'YDir','reverse');
figure(3); plot(1:N,iR,1:N,R);
end
%% save results
save(['HighFPS_GT_',fd.fn], 'V', 'POS', 'R','-v7.3');
end