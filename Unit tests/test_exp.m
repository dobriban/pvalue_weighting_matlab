%Unit tests for Exp Weighting
%%
cd('C:/Git/pvalue_weighting_matlab/Unit tests')
addpath '../Code'

%% negative means
J = 200;
mu = -abs(randn(J,1));
beta = 4;
q = 0.01; %0.05;
w = exp_weights(mu,beta,q);
figure,
plot(mu,w,'*')
filename = 'exp_weights_neg.png';
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
%% any means
J = 200;
mu = randn(J,1);
beta = 4;
q = 0.01; %0.05;
w = exp_weights(mu,beta,q);
figure,
plot(mu,w,'*')
filename = 'exp_weights.png';
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
%% any means - bounded
J = 200;
mu = randn(J,1);
beta = 4;
q = 0.05;
w = exp_weights(mu,beta,q);
figure,
plot(mu,w,'*')
filename = 'exp_weights_bdd.png';
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);