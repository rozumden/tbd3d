function [hmask] = make_hmask(mbb_par, bb, sz)
mbb_pnts = double(mbb_par - [bb(1:2)]'+1);
hmask = logical(poly2mask(mbb_pnts(1,:), mbb_pnts(2,:), sz(1), sz(2)));