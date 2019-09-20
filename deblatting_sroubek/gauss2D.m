function R = gauss2D(mean,sigma,X,Y)
% draw 2D-Gaussian N(mean,sigma) on X,Y mesh
% sigma = [sigma_x, sigma_y] or [sigma] (same for both directions)
% mean = [mean_x, mean_y]
% X ... row vector of x positions of pixel centers
% Y ... column vector of y positions of pixel centers

Xb = (X-0.5-mean(1))/(sqrt(2)*sigma(1));
Xe = (X+0.5-mean(1))/(sqrt(2)*sigma(1));

Yb = (Y-0.5-mean(2))/(sqrt(2)*sigma(end));
Ye = (Y+0.5-mean(2))/(sqrt(2)*sigma(end));

R = (erf(Ye)-erf(Yb))*(erf(Xe)-erf(Xb));
%R = R/sum(R(:));

