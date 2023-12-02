
mu=rand(10,1); cc = cov(rand(100,10)); 
DATA=mvnrnd(mu,cc,500); 
[S,AX,BigAx,H,HAx] = plotmatrix_HJ(DATA,[1 0 0 0.4],':.'); 
