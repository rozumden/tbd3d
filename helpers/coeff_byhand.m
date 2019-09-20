function [coeff] = coeff_byhand()
pnt = ginput();
k = size(pnt,1) - 1;
coeff = {};
for ki = 1:k
	cf = [];
	dr = pnt(ki+1,:) - pnt(ki,:);
	cf = [pnt(ki,:)' dr'];
	coeff = [coeff {cf}];
end