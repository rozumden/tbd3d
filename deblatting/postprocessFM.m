function [F,M] = postprocessFM(F,M)
inl = ~(sum(F,3) < 0.05 | M < 0.05);
inl3 = repmat(inl,[1 1 3]);
M3 = repmat(M, [1 1 3]);
F(inl3) = F(inl3)./M3(inl3);

F(F < 0) = 0;
M(M > 1) = 1;
M(M < 0) = 0;
