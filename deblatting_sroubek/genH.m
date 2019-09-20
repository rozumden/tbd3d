function H = genH(type,param)

% Generates homography matrix 3x3 
% which does any combination of the following depending on "type" 
% type is a cell array of strings or one single string
% and param is a corresponding cell array of parameter vectors
% 'trans' ... translation by param=[tx, ty]
% 'rot' ... rotation by angle param=[a] in degrees
% 'scale' ... scaling param=[sx sy]
% 'H' ... general homography param=[3x3 matrix]


if ischar(type)
   type = {type};
   param = {param};
end

H = eye(3);
for i = 1:length(type)
   switch type{i}
      case 'trans' 
         H = H*[1 0 param{i}(1); 0 1 param{i}(2); 0 0 1];
      case 'rot'
         a = param{i}/180*pi;
         H = H*[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0  0 1];
      case 'scale'
         H = H*[param{i}(1) 0 0; 0 param{i}(2) 0; 0 0 1];
      case 'H'
         H = H*param{i};
   end
end

