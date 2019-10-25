function matF = matF_WB(matF)
fc = [1.9 1 1.8];
WB = [2 1 2]; gamma_coef = 0.4;

for k = 1:3, matF(:,:,k,:) = matF(:,:,k,:) ./ fc(k); end
matF = ((matF.*reshape(WB,1,1,[])/(max(WB))).^gamma_coef);
