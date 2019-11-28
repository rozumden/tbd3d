function [Vk,Vk_WB,PAR,V_WB,V_WB_POS] = generate_lowFPSvideo(V,POS,R,k,resz)
%% input: V,POS,R
% fps/k
if ~exist('k','var')
    k = 8;
end
% downsample low fps video
DSfactor = 1;
% white balance & gamma correction
WB = [2 1 2]; gamma_coef = 0.4;

N = size(V,4);

%% Generate high fps video with WB
V_WB = zeros(size(V),'uint8');
for i = 1:N
    T = double(V(:,:,:,i));
    T = ((T.*reshape(WB,1,1,[])/(255*max(WB))).^gamma_coef)*255;
    V_WB(:,:,:,i) = uint8(T);  
end

%% Generate video with detections
V_WB_POS = zeros(size(V_WB),'uint8');
for i = 1:N
    %figure(1); 
    %dispIm(V_WB(:,:,:,i));
    %hold on; 
    %viscircles(POS(:,i).', R(i),'EdgeColor','b');
    %F = getframe(gca);
    %V_POS(:,:,:,i) = F.cdata;
    if isnan(POS(1,i))
        V_WB_POS(:,:,:,i) = V_WB(:,:,:,i);
    else
        V_WB_POS(:,:,:,i) = insertShape(V_WB(:,:,:,i),'circle',[POS(:,i); R(i)].');
    end
end

%% generate different exposition time k-times longer (Vk)
Vk = zeros(floor(size(V,1)/DSfactor), floor(size(V,2)/DSfactor), 3, floor(N/k),'uint8');
PAR = struct; %zeros(1,floor(N/k));
% fill missing entries in POS by interpolation
inan = find(isnan(POS(1,:)));
fPOS = POS/DSfactor;
if ~isempty(inan)
    i = 1:size(fPOS,2);
    fPOS(1,inan) = interp1(setdiff(i,inan),fPOS(1,setdiff(i,inan)),inan);
    fPOS(2,inan) = interp1(setdiff(i,inan),fPOS(2,setdiff(i,inan)),inan);
end
j = 1;
for i = 1:size(Vk,4)
    % i
    % render image
    e = j+k-1;
    % linear combination without gamma correction
    T = sum(double(V(:,:,:,j:e)),4)/k;
    T = convn(T,ones(DSfactor)/DSfactor^2,'valid');
    Vk(:,:,:,i) = uint8(T(1:DSfactor:end,1:DSfactor:end,:));
    PAR(i).POS = fPOS(:,j:e);
    PAR(i).R = R(j:e);
    j = j+k;
end

Vk_WB = zeros(size(Vk),'uint8');
for i = 1:size(Vk,4)
    T = double(Vk(:,:,:,i));
    T = ((T.*reshape(WB,1,1,[])/(255*max(WB))).^gamma_coef)*255;
    Vk_WB(:,:,:,i) = uint8(T);  
end

if exist('resz','var') && resz ~= 1
    Vk = imresize(Vk, resz);
    Vk_WB = imresize(Vk_WB, resz);
    V_WB = imresize(V_WB, resz);
    for i = 1:numel(PAR)
        PAR(i).R = PAR(i).R * resz;
        PAR(i).POS = PAR(i).POS * resz;
    end
end