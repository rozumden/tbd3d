%% LABS Video Analysis (in a gym)
% Videos acquired in the gym
%
% Position, angle of rotation, axis of rotation and radius are calculated
% from physical properties of a ball rolling on a floor.
%
% Nonlinear lens distortions should be removed, otherwise ball radius is
% not correctly calculated. Videos prior to analysis are therefore rectified using
% Camera Calibibration Toolbox (by Bouguet). - files with name GTundistorted. 
% see ./calibration/Basler2 
%
% Input: load one mat file generated from bounce_test (see code below),
% which contains the video (V) and ball position (POS) and radius
% (R) estimated by Hough transform.
%
% Output: vectors of the video length (values for each frame)
% Omega_GT ... rotation speed (rad/frame)
% LEN_GT ...  length of the ball trajectory (in cm)
% Rx_GT ... radius along x direction (in pixels)
% Ry_GT ... radisu along y durection (in pixels)
% R_GT ... mean of Rx and Ry (in pixels)
% p{1} ... ball center (in pixels)
% p{2} ... point halfway between the point of contact with floor and
% rotation axis (in pixels); good for visually checking in the video
% that the ball in exactly in the same position after rotating 360degrees. 
% p{3} ... position of the rotation axis leaving the ball surface (in pixels)
% p{4} ... left, right, top, and bottom extremal point on the ball; good
% for checking accuracy of the estimated radius
% 
% Note: For validating our TbD approach, we should be comparing with Omega_GT, p{3}
% and R_GT

%% lab1
fname = 'HighFPS_GT_lab2';
load(fname);
% trajectory parameters
Lp = [89; 215]; %[87; 213]; %Lpoint in the image [pixels]
Rp = [807; 205]; %[807; 203]; % Rpoint in the image [pixels]
Ldist = [317.5;-85.8]-[0;0]; %cm
Rdist = [313.5;73]-[0;0]; %cm
Cdist = 305.1 + 5; %cm
%% lab2_Rectified
load HighFPS_GTundistorted_lab2.mat
Lp = [84; 213];
Rp = [809; 203];
Ldist = [317.5;-85.8]-[5;0]; %cm
Rdist = [313.5;73]-[5;0]; %cm
Cdist = 300.1 + 5; %305.5; %cm %optical axis is slighlty shifted from the center; see Basler.cc
%% lab3d_1 Rectified
load HighFPS_GTundistorted_lab3D_1.mat
Lp = [79; 289];
Rp = [671; 177];
Ldist = [245.5; -85.8]-[5; 0];
Rdist = [348; 73]-[5;0];
Cdist = 295+7; %294; %optical axis is slighlty shifted from the center; see Basler.cc
%% lab3d_1
load HighFPS_GT_lab3D_1.mat
Lp = [83; 288];
Rp = [666; 177];
Ldist = [245.5; -85.8]-[0;0];
Rdist = [348; 73]-[0;0];
Cdist = 306.9+7; %294; %optical axis is slighlty shifted from the center; see Basler.cc
%% lab3d_1 Rectified2
fname = 'HighFPS_GTundistorted2_lab3D_1'
load(fname);
Lp = [74; 290];
Rp = [669; 177];
Ldist = [245.5; -85.8]-[9; 0];
Rdist = [348; 73]-[9;0];
Cdist = 294.7+7; %294; %optical axis is slighlty shifted from the center; see Basler.cc
%% lab2_Rectified2
fname = 'HighFPS_GTundistorted2_lab2'
load(fname);
Lp = [82; 214];
Rp = [810; 204];
Ldist = [317.5;-85.8]-[4.5;0]; %cm
Rdist = [313.5;73]-[4.5;0]; %cm
Cdist = 300.5 + 5; %305.5; %cm %optical axis is slighlty shifted from the center; see Basler.cc


%% 
% ball parameters
ball_circum = 74.9; %cm
ball_radius = ball_circum/(2*pi);
% ball was guided by side tracks
track_width = 10; %cm
h = ball_radius-track_width/2;
cap_radius = sqrt(h*(2*ball_radius-h));
cap_circum = 2*pi*cap_radius;

% camera parameters
camera_height = 55.5; %cm
%camera_f = 1.6; %cm
Basler = load('./calibrate/Basler/Calib_Results');
% camera intrinsic matrix
KK = [Basler.fc(1) Basler.alpha_c*Basler.fc(1) Basler.cc(1);0 Basler.fc(2) Basler.cc(2) ; 0 0 1];
% correction: cropped (from 200 to 1000 line), raw (half size), and left top pixel is [1 1]
T =  [1 0 1; 0 1 1; 0 0 1]*[0.5 0 0; 0 0.5 0; 0 0 1]*[1 0 0; 0 1 -200; 0 0 1];


L3 = (T*KK)\[Lp; 1];
%L(1:2) = apply_distortion(L(1:2),Basler.kc);
R3 = (T*KK)\[Rp; 1];
%R(1:2) = apply_distortion(R(1:2),Basler.kc);
L3 = L3*norm([Ldist(1),camera_height])/norm(L3);
R3 = R3*norm([Rdist(1),camera_height])/norm(R3);
disp(['Distance (L,R) differs from the measured by ',num2str(norm(L3-R3)-(Rdist(2)-Ldist(2))),'cm.']);
%nv = -cross(L3,R3)/norm(cross(L3,R3));
tv = (R3-L3)/norm(R3-L3);
nv = cross(tv,[0;0;1])/norm(cross(tv,[0;0;1]));
phi = atan(camera_height/Cdist);
RM = cos(phi)*eye(3) + sin(phi)*[0 -tv(3) tv(2); tv(3) 0 -tv(1); -tv(2) tv(1) 0] ...
    + (1-cos(phi))*(tv*tv.');
shift = ball_radius*RM*nv;
% shift L and R point to be in the center of the ball
L3 = L3+shift;
R3 = R3+shift;
Lc = T*KK*L3; Lc = Lc(1:2)/Lc(3);
Rc = T*KK*R3; Rc = Rc(1:2)/Rc(3);

v = Rc-Lc;
s = v.'*(POS-Lc)/(v.'*v);
t = s*L3(3)./(s*(L3(3)-R3(3))+R3(3));
y = t*sqrt(sum((R3-L3).^2,1));
% fit parabolic curve (assumption of constant decceleration)
y = y(:);
x = (1:length(y)).';
A = [x.^2, x, ones(size(x))];
a = (A.'*A)\(A.'*y);
LEN_GT = (a(1)*x.^2 + a(2)*x + a(3)).';
% correct back 
t = LEN_GT/sqrt(sum((R3-L3).^2,1));
s = t*R3(3)./(L3(3)-t*(L3(3)-R3(3)));
d = L3(3)./(s*(L3(3)/R3(3) - 1) + 1);
% projected ball center point
sp = (R3-L3)*t + L3;
p{1} = v.*s + Lc;
dt = 1;
Omega_GT = 360*(LEN_GT(1+dt:end)-LEN_GT(1:end-dt))/cap_circum/dt; % in degree
Omega_GT = [Omega_GT, Omega_GT(end)];
% projected ball point half way from front to bottom
phi2 = pi/4;
RM2 = cos(phi2)*eye(3) + sin(phi2)*[0 -tv(3) tv(2); tv(3) 0 -tv(1); -tv(2) tv(1) 0] ...
    + (1-cos(phi2))*(tv*tv.');
shift = RM2*(ball_radius*cross(tv,RM*nv)); %-ball_radius*nv;
p{2} = T*KK*(sp + shift);
p{2} = p{2}(1:2,:)./p{2}(3,:);
% projected ball front point (rotation axis)
shift = ball_radius*cross(tv,RM*nv);
p{3} = T*KK*(sp + shift);
p{3} = p{3}(1:2,:)./p{3}(3,:);
% projected ball bottom point (point in contact with the floor)
%shift = -ball_radius*RM*nv;
%p4 = T*KK*((R-L).*t + L + shift);
%p4 = p4(1:2,:)./p4(3,:);

%sp = (T*KK)\[POS; ones(1,size(POS,2))];
%sp = sp.*d;
n = cross(repmat([0; 1; 0],1,size(t,2)),sp);
n = n./sqrt(sum(n.^2,1));
shift = ball_radius*n;
p{4}.R = T*KK*(sp + shift);
p{4}.R = p{4}.R(1:2,:)./p{4}.R(3,:);
p{4}.L = T*KK*(sp - shift);
p{4}.L = p{4}.L(1:2,:)./p{4}.L(3,:);
n = cross(repmat([1; 0; 0],1,size(t,2)),sp);
n = n./sqrt(sum(n.^2,1));
shift = ball_radius*n;
p{4}.T = T*KK*(sp + shift);
p{4}.T = p{4}.T(1:2,:)./p{4}.T(3,:);
p{4}.B = T*KK*(sp - shift);
p{4}.B = p{4}.B(1:2,:)./p{4}.B(3,:);

% due to the camera perspective, the ball is deformed: different radius along x and y
% we take the mean (is this correct?)
Rx_GT = sqrt(sum((p{4}.L-p{4}.R).^2,1))/2;
Ry_GT = sqrt(sum((p{4}.T-p{4}.B).^2,1))/2;
R_GT = mean([Rx_GT; Ry_GT]);

%%
figure; plot(1:length(R),R,1:length(R),R_GT,1:length(R),Rx_GT,1:length(R),Ry_GT);
%%
figure; dispIm(linscale(double(V_WB(:,:,:,1))));
hold on;
plot([Lp(1);Rp(1)],[Lp(2);Rp(2)]);

%% Save results
save([fname,'_Phys'],'-v7.3','Omega_GT','p','LEN_GT','Rx_GT','Ry_GT','R_GT')
%% visually check the detection
VPOS = V_WB;
for i = 1:N
    if prod(POS(:,i)==0)
        continue;
    end
    %VPOS(POS(2,i),POS(1,i),:,i) = repmat(cat(3,255,0,0),1,1);
    %VPOS(round(SPOS(2,i)),round(SPOS(1,i)),:,i) = repmat(cat(3,0,255,0),1,1);
    VPOS(round(p{1}(2,i)),round(p{1}(1,i)),:,i) = repmat(cat(3,255,0,0),1,1);
    VPOS(round(p{2}(2,i)),round(p{2}(1,i)),:,i) = repmat(cat(3,0,255,0),1,1);
    VPOS(round(p{3}(2,i)),round(p{3}(1,i)),:,i) = repmat(cat(3,255,255,0),1,1);
    if sum(round(p{4}.L(:,i))>1) == 2 
        VPOS(round(p{4}.L(2,i)),round(p{4}.L(1,i)),:,i) = repmat(cat(3,0,0,255),1,1);
        VPOS(round(p{4}.R(2,i)),round(p{4}.R(1,i)),:,i) = repmat(cat(3,0,0,255),1,1);
        VPOS(round(p{4}.T(2,i)),round(p{4}.T(1,i)),:,i) = repmat(cat(3,0,0,255),1,1);
        VPOS(round(p{4}.B(2,i)),round(p{4}.B(1,i)),:,i) = repmat(cat(3,0,0,255),1,1);
    end
end
%% Generate high fps video with WB
% white balance & gamma correction
WB = [2 1 2]; gamma_coef = 0.4;
N = size(V,4);
V_WB = zeros(size(V),'uint8');
for i = 1:N
    T = double(V(:,:,:,i));
    T = ((T.*reshape(WB,1,1,[])/(255*max(WB))).^gamma_coef)*255;
    V_WB(:,:,:,i) = uint8(T);  
end
%%

for i = 1:size(SPOS,2)
    if mod(i,1)==0
        figure(1); clf;
        dispIm(linscale(double(V_WB(:,:,:,i))));
        hold on;
        p = v*s(i) + Lp;
        plot([SPOS(1,i); p(1)],[SPOS(2,i); p(2)],'x');
        keyboard;
    end
end

%% find the closest point on the path to the center point
L1 = (T*KK)\[Lp; 1];
R1 = (T*KK)\[Rp; 1];
L1 = L1*norm([Ldist(1),camera_height])/norm(L1);
R1 = R1*norm([Rdist(1),camera_height])/norm(R1);
u1 = R1-L1;
b1 = L1;
v1 = [0; 0; 1];
s1 = (-u1'*b1 + ((u1'*v1)*(v1'*b1))/(v1'*v1))/(-((u1'*v1)*(v1'*u1))/(v1'*v1)+u1'*u1);
C1 = u1*s1+b1;
sqrt(C1(3)^2-camera_height^2)
a1 = (T*KK)*v1; a1 = a1(1:2)/a1(3);
b1 = (T*KK)*C1; b1 = b1(1:2)/b1(3);
%% quadratic regression
y = LEN_GT;
y = y(:);
x = (1:length(y)).';
A = [x.^2, x, ones(size(x))];
a = (A.'*A)\(A.'*y);
y = a(1)*x.^2 + a(2)*x + a(3);
%% linear regression to P
x = (1:length(L3)).';
A = [x,ones(length(L3),1)];
a = (A'*A)\(A'*L3.');
y = a(1)*x + a(2);
figure; plot(x,L3,x,y);