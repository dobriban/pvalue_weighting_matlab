%Unit tests for Exp Weighting
%%
cd('C:/Git/pvalue_weighting_matlab/Unit tests')
addpath '../Code'

%%
J = 200;
mu = -abs(randn(J,1));
beta = 4;
q = 0.01; %0.05;
w = exp_weights(mu,beta,q);
plot(mu,w,'*')

