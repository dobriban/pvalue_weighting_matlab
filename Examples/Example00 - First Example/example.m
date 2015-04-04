%Example for the P-value Weighting Code for Matlab 
%%
addpath '../../Code'
addpath '../../Code/Helper Code'
addpath '../../Data'

%% Step I: Generate Data or load data
%% Option (1): load readily generated data
load('../../Data/example_data.mat', 'J','P_current','t1','t2');

%% Option (2): generate data
%generate data from a mixture distribution
%mu_i ~ -absolute N(0,1)
%X_i ~ Bernoulli(p)
%(t1_i,t2_i) ~ Normal(X(mu_i,mu_i), Id_2)
%P_prior(i) = Phi

rng(0); %set seed
J = 1e3;
mu = - 2*abs(randn(J,1));
frac_sig = 0.1;
X = binornd(1,frac_sig,J,1);
t1 = normrnd(X.*mu,1);
t2 = normrnd(X.*mu,1);
P_current = normcdf(t2); %current p-values
%%
save('../../Data/example_data.mat', 'J','P_current','t1','t2');

%% Visualize
scatter(t1,t2);
xlabel('Prior Test Statistic');
ylabel('Current Test Statistic');
title('Synthetic Data Example');
saveTightFigure(gcf,'../../Data/Synthetic_Data_example.pdf')


%% data analysis
%% (1) unweighted
alpha = 0.05; report = 'yes';
h_u=bonferroni(P_current,alpha,report);

%% (2) Regularized weighting
q =alpha/J; %expected fraction of false rejections under 'null'
sigma = ones(J,1);
w = bayes_weights(t1,sigma,q);
P_wr = P_current./w;
h_r=bonferroni(P_wr,alpha,report);

%% Post-Analysis
find(h_u==1)
find(h_r==1)

%%
plot(t1,w,'*');
xlabel('Prior Test Statistic');
ylabel('weight');
title('Plot of Weights');
saveTightFigure(gcf,'../../Data/weights.pdf')

