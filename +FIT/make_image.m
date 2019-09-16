function [Hvis] = make_image(Hv, curve, curve_full, st, en)
alp = 0.5;
alp_endpoints = 0;
clr_full = [0 1 0];
clr = [1 0 0];
clr_st = [0 1 0];
clr_en = [0 0 1];
N = size(Hv, 1) * size(Hv, 2);
Hvis = repmat(Hv, [1 1 3/size(Hv,3)]);
Hvis0 = Hvis;
for k = 1:3
	if nargin > 2
		Hvis((k-1)*N+sub2ind(size(Hv), curve_full(2,:), curve_full(1,:))) = ...
				alp*Hvis(sub2ind(size(Hv), curve_full(2,:), curve_full(1,:))) + (1-alp)*clr_full(k);  
	end
	
	Hvis((k-1)*N+sub2ind(size(Hv), curve(2,:), curve(1,:))) = ...
			alp*Hvis0(sub2ind(size(Hv), curve(2,:), curve(1,:))) + (1-alp)*clr(k); 

	if nargin > 3
		st = round(st);
		Hvis(st(2),st(1),k) = alp_endpoints*Hvis(st(2),st(1),k) + (1-alp_endpoints)*clr_st(k);  
	end 
	if nargin > 4
		en = round(en);
		Hvis(en(2),en(1),k) = alp_endpoints*Hvis(en(2),en(1),k) + (1-alp_endpoints)*clr_en(k);  
	end 
end 
