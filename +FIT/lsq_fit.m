function [coeff] = lsq_fit(frms, st_fit, en_fit, max_power, expos)
if nargin < 5
	expos = 1;
end
if nargin < 4
	max_power = 2;
end
A = [];
b = [];
Aeq = [];
beq = [];
for ifit = st_fit:en_fit
	frm = frms{ifit};
	if isempty(frm), continue; end	
	if numel(frm.Start) ~= 2 || numel(frm.End) ~= 2, error('No Start or End'); end
	Ast = [];
	Aen = [];
	for pow1 = 0:max_power
		Ast=[Ast [ifit.^pow1 0; 0 ifit.^pow1]];
		Aen=[Aen [(ifit+expos).^pow1 0; 0 (ifit+expos).^pow1]];
 	end
 	
	if ifit == st_fit
		Aeq = [Aeq; Ast];
		beq = [beq; frm.Start'];
	else
		A = [A; Ast];
		b = [b; frm.Start'];
	end

	if ifit == en_fit
		Aeq = [Aeq; Aen];
		beq = [beq; frm.End'];
	else
		A = [A; Aen];
		b = [b; frm.End'];
	end
end
if isempty(A)
	A = Aeq;
	b = beq;
end

if isempty(A)
	coeff = {ones(2,2)};
	return;
end

xs = A \ b;
nvars = (max_power+1)*2;
options = optimset('Display','none');
xs = lsqlin(A,b,[],[],Aeq,beq,-Inf*ones(nvars,1),Inf*ones(nvars,1),[],options);

coeff = {reshape(xs, [2 max_power+1])};
