function X = from2Dto3D(x,r)

Basler = load('Calib_Results');
% camera intrinsic matrix
KK = [Basler.fc(1) Basler.alpha_c*Basler.fc(1) Basler.cc(1);0 Basler.fc(2) Basler.cc(2) ; 0 0 1];
% correction: cropped (from 200 to 1000 line), raw (half size), and left top pixel is [1 1]
T =  [1 0 1; 0 1 1; 0 0 1]*[0.5 0 0; 0 0.5 0; 0 0 1]*[1 0 0; 0 1 -200; 0 0 1];

P = (T*KK)\[x(:); 1];

X = P/norm(P)*(1/r);