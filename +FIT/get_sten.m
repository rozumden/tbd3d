function [st,en] = get_sten(frm)
coeff0 = cellfun(@(x)fliplr(x).', frm.coeff, 'UniformOutput', false);
for c = 1:numel(coeff0)
	if ~isempty(frm.bb)
		coeff0{c}(:,1) = coeff0{c}(:,1) + [frm.bb(1:2)]' - 1;
	end
end

pnts = FIT.trajPredictor(coeff0, frm.Direction, 1)';
en = pnts(:,1);

pnts = FIT.trajPredictor(coeff0, -frm.Direction, 1)';
st = pnts(:,1);



