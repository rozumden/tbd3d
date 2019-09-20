function [H, F, M, report] = convsparseCG_v1z(G,B,F,M,H,par)
% Solve the deblatting problem
% Find object trajetories, appearance and shape from a single image and
% background
%
% [H, F, M, report] = convsparseCG_v2(G,B,F,M,H,par)
%
% input:
% G ... color blurred image MxNxL (L... number of color channels)
% B ... color background image MxNxL (L... number of color channels)
% F ... filters (color foreground objects); cell array { (M1xN1xL) (M2xN2xL) ...}
% M ... masks; cell array in the same format as F; default is { 0 0 ...}
%       if cells contain scalars 0 or 1: 0... generate masks full of ones; 1... get masks from F (F~=0)  
% H ... initial PSFs (trajectories of objects F) MxNxlength(F) (default is
%       empty)
% par . structure containing parameters that override default ones
%
% output:
% H ... estimated PSFs (MxNxlength(F))
% F ... estimated foreground objects (object appearance)
% M ... estimated masks (object shape)
% report ... structure; 
%            field E: energy in every iteration; 
%            field PAR: parameters of the algorithm
%
% Note: Same as convsparseCG_v1 but with TV regularization also across
% individual objects in F and M (along the third dimension). Usefull if we
% expect small changes in F and M along the trajectory. 
% Note: PrM = asetupLnormPrior(Lp,0.1*betaF,betaF); in _v1 it is 0.5*betaF
% Note: minHstep is uncommented 

%% DEFAULT PARAMETERS
% Number of iterations in each step
maxiter_main = 100; % max number of iterations in the main loop
maxiter_fm = 0; % number of initial minF and minM iterations
maxiter_H = 1; % max number of iterations in the minHstep
maxiter_F = 1; % max number of iterations in the minFstep
maxiter_M = 1; % max number of iterations in the minMstep

VERBOSE = 2;
%% MAIN PARAMETERS
%(alpha, beta, xi, maxbeta, tau, gamma)
% 1e1 1e1 1e1 1e3 1.1 1e0
% 1e1 1e2 1e2 1e4 1.1 1e1
Lp = 1; % type of sparse norm l_p
alpha = 1e1; %1; % sparsity weight for H
alphaF = 1e-3; % 1e-2; % 1e-4; % TV(F) weight 
beta = 1e3; %1e2; %1e3; %1e1; % beta || a - v ||^2
betaF = 1e-1; %1e-1; %1e0; %1e-2; 
betaM = 1e-1; %1e4; %1e2; % weight for projection to FM simplex
%maxbeta = 1e1; %1e3;
%tau = 1.1;  
gamma = 1e2; %1e0 for pinkpong; %1e-4 was for image intensities 0..255 % data term weight

reltol = 1e-4; % CG reltol
%% DEFAULT PARAMETERS END HERE

if exist('par','var') && ischar(par)
    [pathstr,name] = fileparts(par);
    curf = pwd;
    if ~isempty(pathstr)
        cd(pathstr);
    end
    eval(name);
    cd(curf);
end

if nargin < 4 || isempty(M)
    flag = num2cell(zeros(1,length(F)));
else
    flag = M;
end

gsize = [size(G,1),size(G,2),size(G,3)];
L = gsize(3); % number of color channels

if VERBOSE>1
    FIG_HANDLE_F = figure(1);
    FIG_HANDLE_H = figure(2);
    %FIG_HANDLE_M = figure(3);
end


% variables used in min_H
Hsize = [gsize(1),gsize(2),length(F)];
hsize = [gsize(1)*gsize(2),length(F)];
Tsize = [Hsize, L];
tsize = [hsize, L];
Fsize = [size(F{1},1), size(F{1},2), length(F), L];
Msize = Fsize(1:3);
HTH = 0;
% shift H as we use fft2 instead of conv to perform convolution
H = circshift(H,-floor(Msize(1:2)/2));
if nargin < 5 || isempty(H)
    H = zeros(Hsize);
    FH = zeros(hsize);
else
    FH = reshape(fft2(H),hsize);
end
V = zeros(Hsize); % Auxiliar variable
A = zeros(Hsize); % Lagrange multiplier
Pr = asetupLnormPrior(Lp,alpha,beta);
V = Pr.fh(H-A,abs(H-A));

% variables used in min_F
Dx = [1 -1; 0 0];
Dy = [1 0; -1 0];
if Fsize(3) > 1
    Dz = cat(3,1,-1);
else
    Dz = []; % if only one F then do not use TV across different F (z axis)
end
%Dz = [];
DTD2 = conv2(Dx,rot90(Dx,2),'full') + conv2(Dy,rot90(Dy,2),'full');
W = zeros(Fsize); % Auxiliar variable
C = zeros(Fsize); % Lagrange multiplier
Dsize = Fsize; 
Dsize(1:2) = Dsize(1:2) + 1;
Wx = zeros(Dsize);
Cx = zeros(Dsize);
Wy = zeros(Dsize);
Cy = zeros(Dsize);
Wz = zeros(Dsize);
Cz = zeros(Dsize);

D = zeros(Msize); % Lagrange multiplier
Z = zeros(Msize); % Auxiliar variable

Zx = zeros(Dsize(1:3));
DDx = zeros(Dsize(1:3));
Zy = zeros(Dsize(1:3));
DDy = zeros(Dsize(1:3));
Zz = zeros(Dsize(1:3));
DDz = zeros(Dsize(1:3));

if length(flag) == 1 % one mask for all F's
    disp('Using one mask.');
    one_mask = true;
    flag = repmat(flag,1,Fsize(3));
else
    one_mask = false;
end
if length(flag) ~= Fsize(3)
    error('Size of F and M does not match!');
end
for c = 1:hsize(2)
    %FF(:,c,:) = reshape(fft2(double(F{c}),gsize(1),gsize(2)),[hsize(1), 1, L]);
    if numel(flag{c}) ==  1      
        if flag{c} == 1 %initialize mask with ones or as delta
            M{c} = double(any(F{c}~=0,3));
        else
            M{c} = ones(size(F{c},1),size(F{c},2));
        end
    else %initialize mask with the input argument
        M{c} = flag{c}; 
    end
    %FM(:,c,:) = reshape(fft2(double(M{c}),gsize(1),gsize(2)),[hsize(1), 1, L]);
end
M = reshape(cell2mat(M),Msize); % replicate mask for every color channel
F = reshape(cell2mat(F),Fsize);

% set F to delta
%F = zeros(Fsize);
%F(floor(Fsize(1)/2)+1, floor(Fsize(2)/2)+1,:, :) = 1;

% set F to ones
%F = zeros(Fsize); F(:,:,:,3) = 1;

%PrF = asetupLnormPrior(Lp,alphaF,betaF);
%W = PrF.fh(F-C,abs(F-C));
W = F; W(W<0) = 0;
Z = M; Z(Z<0) = 0; Z(Z>1) = 1;
FF = reshape(fft2(F,gsize(1),gsize(2)),tsize);
FM = reshape(fft2(M,gsize(1),gsize(2)),hsize);

% common input variables
FGB = reshape(fft2(G-B),[hsize(1),1,L]);
%rB = B(Fsize(1):end,Fsize(2):end,:);
%rGB = G(Fsize(1):end,Fsize(2):end,:)-rB;


report.E = [];

fh_projection = @(x,y) (project2pyramid(x,y,1)); % projection to a convex set of arbitrary dimension
% initial estimation of F
for j = 1:maxiter_fm
    minFMstep;
end

report.PAR = struct('Lp', Lp, 'alpha', alpha, 'alphaF', alphaF, 'beta', beta, ...
    'betaF', betaF, 'betaM', betaM, 'gamma', gamma, 'reltol', reltol);
report.E = [report.E; calcE];

for i_main = 1:maxiter_main
    if VERBOSE
        disp([num2str(i_main)]);
    end
    minFMstep;
    
    %minHstep;
    
    report.E = [report.E; calcE];
end

%shift H back for use with conv instead of fft2
H = circshift(H,floor(Msize(1:2)/2));



%% END

%% update H (object trajectory)
function minHstep

Fb1 = sum(conj(FF).*FGB,3) - ...
        sum(conj(FM).*reshape(fft2(B.*(G-B)),[hsize(1),1,L]),3);
b1 = real(ifft2(reshape(Fb1,Hsize)));
% additive
Pr = asetupLnormPrior(Lp,alpha,beta); 

for i=1:maxiter_H
    % size hsize 
    b = b1 + beta/gamma*(V+A);
    
    [xmin,fl,relres,iter,resvec] = mycg(@gradH,vec(b),reltol,1000,[],vec(H));
    if VERBOSE
        disp(['minH ',num2str(i),' - flag, iter:',num2str([fl iter])]);
    end
    H = reshape(xmin,Hsize);
    
    
    %% additive
    V = Pr.fh(H-A,abs(H-A)); % Lp norm 
    V(V<0) = 0;
    % update Bregman variable (Lagrange variable)
    A = A + V - H;
    if VERBOSE>1
        update4Fig(FIG_HANDLE_H,{'V', reshape(V,Hsize(1),[])},{'H', reshape(H,Hsize(1),[])},{'A', reshape(A,Hsize(1),[])},...
            {'G-B',linscale(G-B)} );   
    end
end
FH = reshape(fft2(H),hsize);
%report.E = [report.E, calcE];
end

function r = gradH(x)
    %X = reshape(x,hsize);
    X = reshape(fft2(reshape(x,Hsize)),hsize);
    % (F - diag(B) M)
    FT = reshape(sum(FF.*X,2),[gsize(1),gsize(2),L]) - ...
        fft2(B.*real(ifft2(reshape(sum(FM.*X,2),[gsize(1),gsize(2)]))));
    T = real(ifft2(FT));
    % (F^T - M^T diag(B))
    FT = reshape(FT,[hsize(1) 1 L]);
    T = sum(conj(FF).*FT,3) - ...
        sum(conj(FM).*reshape(fft2(B.*T),[hsize(1),1,L]),3);
    % beta/gamma I
    %r = vec(T) + beta/gamma*x;
    r = vec(real(ifft2(reshape(T,Hsize)))) + beta/gamma*x;
end

%% update F and M simultenously
function minFMstep
b1 = conj(FH).*FGB;
b1 = real(ifft2(reshape(b1,[gsize(1), gsize(2), hsize(2), L])));
b1 = b1(1:Fsize(1),1:Fsize(2),:,:);
b2 = conj(FH).*sum(reshape(fft2(B.*(G-B)),size(FGB)),3);
b2 = real(ifft2(reshape(b2,[gsize(1), gsize(2), hsize(2)])));
b2 = -b2(1:Fsize(1),1:Fsize(2),:);
Pr = asetupLnormPrior(Lp,alphaF,betaF); 
PrM = asetupLnormPrior(Lp,0.1*betaF,betaF);
for i=1:maxiter_F
    ba = b1 + betaM/gamma*(W+C) + ...
        betaF/gamma*(convn(Wx+Cx,rot90(Dx,2),'valid') + convn(Wy+Cy,rot90(Dy,2),'valid') + zconvn(Wz+Cz,Dz,'t'));
    bb = b2 + betaM/gamma*(Z+D) + ...
        betaF/gamma*(convn(Zx+DDx,rot90(Dx,2),'valid') + convn(Zy+DDy,rot90(Dy,2),'valid') + zconvn(Zy+DDy,Dz,'t'));
    [xmin,fl,relres,iter,resvec] = mycg(@gradFM,[ba(:);bb(:)],reltol,1000,[],[F(:);M(:)]);
    
    F = reshape(xmin(1:prod(Fsize)),Fsize);
    M = reshape(xmin(prod(Fsize)+1:end),Msize);
    
    W = F-C;
    if one_mask
        T1 = mean(M,3) - D(:,:,1);
        [T1, W(:)] = fh_projection(T1(:), reshape(W,[],Fsize(3)*L));
        Z(:) = repmat(T1,Fsize(3),1);
    else
        Z = M-D;
        [Z(:), W(:)] = fh_projection(Z(:), reshape(W,[],L));
    end
   
       
    xD = convn(F,Dx,'full');
    yD = convn(F,Dy,'full');
    zD = zconvn(F,Dz,'n');
    xDm = xD - Cx;
    yDm = yD - Cy;
    zDm = zD - Cz;
    nDm = repmat(sqrt(sum(xDm.^2,4) + sum(yDm.^2,4) + sum(zDm.^2,4)),[1 1 1 Fsize(4)]);
    Wy = Pr.fh(yDm,nDm);
    Wx = Pr.fh(xDm,nDm);
    Wz = Pr.fh(zDm,nDm);
 
    xM = convn(M,Dx,'full');
    yM = convn(M,Dy,'full');
    zM = zconvn(M,Dz,'n');
    if one_mask
        xM = repmat(mean(xM,3),1,1,size(xM,3));
        yM = repmat(mean(yM,3),1,1,size(yM,3));
        zM = repmat(mean(zM,3),1,1,size(zM,3));
    end
    xDm = xM - DDx;
    yDm = yM - DDy;
    zDm = zM - DDz;
    nDm = sqrt(xDm.^2 + yDm.^2 + zDm.^2);
    Zy = PrM.fh(yDm,nDm);
    Zx = PrM.fh(xDm,nDm);
    Zz = PrM.fh(zDm,nDm);
    
    % update Bregman variables (Lagrange variables)
    D = D - M + Z;
    C = C + W - F;
    Cx = Cx + Wx - xD;
    Cy = Cy + Wy - yD;
    Cz = Cz + Wz - zD;
    
    DDx = DDx + Zx - xM;
    DDy = DDy + Zy - yM;
    DDz = DDz + Zz - zM;
    
    if VERBOSE
        disp(['minFM ',num2str(i),' - flag, iter: ',num2str([fl iter])]);
    end
    if VERBOSE>1
        update4Fig(FIG_HANDLE_F,{'Z', linscale(reshape(Z,Fsize(1),[]))},...
            {'F', linscale(reshape(F,Fsize(1),[],L))},...
            {'M', linscale(reshape(M,Fsize(1),[]))},...
            {'W', linscale(reshape(W,size(W,1),[],L))});     
    end
end
FF = reshape(fft2(F,gsize(1),gsize(2)),tsize);
FM = reshape(fft2(M,gsize(1),gsize(2)),hsize);
end
function r = gradFM(x)
    Xf = reshape(x(1:prod(Fsize)),Fsize);
    Xm = reshape(x(prod(Fsize)+1:end),Msize);
    FXf = reshape(fft2(Xf,Hsize(1),Hsize(2)),tsize);
    FXm = reshape(fft2(Xm,Hsize(1),Hsize(2)),tsize(1:2));
    T = sum(FH.*FXf,2) - reshape(fft2(B.*real(ifft2(reshape(sum(FH.*FXm,2),Hsize(1:2))))),size(FGB));
    Tf = real(ifft2(reshape(conj(FH).*T,Tsize)));
    Tf = Tf(1:Fsize(1),1:Fsize(2),:,:);
    Tm = conj(FH).*reshape(fft2(sum(B.*real(ifft2(reshape(T,gsize))),3)),[hsize(1),1]);
    Tm = real(ifft2(reshape(Tm,[gsize(1), gsize(2), hsize(2)])));
    Tm = -Tm(1:Fsize(1),1:Fsize(2),:);
    r = [vec(Tf + betaM/gamma*Xf + betaF/gamma*(convn(Xf,DTD2,'same') + zconvn(Xf,Dz,'tn')));...
        vec(Tm + betaM/gamma*Xm + betaF/gamma*(convn(Xm,DTD2,'same') + zconvn(Xm,Dz,'tn')))];
end

function E = calcE
    FE = reshape(sum(FF.*repmat(FH,1,1,L),2),[gsize(1),gsize(2),L]) - ...
           fft2(B.*real(ifft2(reshape(sum(FM.*repmat(FH,1,1,L),2),[gsize(1),gsize(2),L])))) - ...
           reshape(FGB,[gsize(1),gsize(2),L]);
    E = sum(abs(FE(:)).^2)/numel(FE);
end

function R = zconvn(A,B,type)
    % function is correct only for B = cat(3,1,-1)
    if isempty(B) % do nothing if we do not want to calculate TV across F and M
        R = 0;
        return;
    end
    switch type
        case 'n'
            R = cat(3,padarray(convn(A,B,'valid'),[1 1],0,'post'), ...
                zeros(size(A,1)+1,size(A,2)+1,1,size(A,4)));
        case 't'
            R = cat(3,-A(:,:,1,:),convn(A(:,:,1:end-1,:),flip(B,3),'valid'),A(:,:,end-1,:));
            R = R(1:end-1,1:end-1,:,:);
        case 'tn'
            BTB = convn(B,flip(B,3),'full');
            R = cat(3, A(:,:,1,:)-A(:,:,2,:), convn(A,BTB,'valid'), -A(:,:,end-1,:)+A(:,:,end,:));
    end
end

function r = flip3(x)
    r = flip(flip(flip(x,1),2),3);
end

function R = myconvn(A,B,shape)
    R = convn(A,B,'full');
    s = [(size(B,1)-1)*shape(1), (size(B,2)-1)*shape(2), (size(B,3)-1)*shape(3), (size(B,4)-1)*shape(4)];
    R = R(1+s(1):end-s(1),1+s(2):end-s(2),1+s(3):end-s(3),1+s(4):end-s(4));
end

end

