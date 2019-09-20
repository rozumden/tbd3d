function [I,curve,X,Y] = genstreak(Rz,T,rho,y0,imgsize)
% generate general streak of point yo given by rotation (Rz) along z in 
% time, translation (T) x,y,z in time and depth (rho)
%
% Rz = [angle_0, angle_1, ....]
% T = [  Tx_0, Tx_1 ...;
%        Ty_0, Ty_1 ...;
%        Tz_0, Tz_1 ...]
% y0 = [x,y] ... point coordinates


% sigma of the Gaussian
sigma = 0.5;

% number of points on the streak
N = size(T,2);
y0 = y0(:);

% determine the streak size
curve = zeros(2,N);


for i = 1:N
  c = genH('rot',Rz(i))*([y0; rho] + T(:,i));
  curve(:,i) = c(1:2)/c(3);  
end


roi = [floor(min(curve.',[],1))-1; ceil(max(curve.',[],1))+1];
I = zeros(imgsize);
X = -floor((size(I,2)-1)/2):ceil((size(I,2)-1)/2);
Y = [-floor((size(I,1)-1)/2):ceil((size(I,1)-1)/2)].';
 
for i = 1:N
   I = I + gauss2D(curve(:,i),sigma,X,Y);
end
I = I/sum(I(:));

end
