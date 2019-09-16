function [curve] = fit_parabola(c0, c1, st, en)
%%    [a b c d e f]
Aeq = [1 0 0 0 0 0; % start point x
	   0 0 0 1 0 0; % start point y
	   1 1 1 0 0 0; % end point x
	   0 0 0 1 1 1; % end point y
   		];
C = [0 1 0 0 0 0; % start derivate x
	 0 0 0 0 1 0; % start derivate y
	 0 1 2 0 0 0; % end derivate x
	 0 0 0 0 1 2; % end derivate y
 ];
beq = double([st; en]);
d = double([ [c0{end}(:,2)+2*c0{end}(:,3)]; c1{1}(:,2)]);
options.Display = 'off';
x = lsqlin(C,d,[],[],Aeq,beq,[],[],[],options);
curve = {[x([1 2 3; 4 5 6])]};