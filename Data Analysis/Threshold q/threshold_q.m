%compute the thresholds for values of q for the three data examples
%see: 0.2447, 0.2836, 0.2337
cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Analysis\Cardiogram-C4D');
addpath '../../Code'
addpath '../../Code/Weighting'
%% Data set 1:
prior = 'Cardiogram-CAD.Nat.Genet.2011';
current = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';
[~, ~, N_prior, Z_prior,~, N_current] = load_data(prior, current);
%%
mu = sqrt(N_current./N_prior).*Z_prior;
sigma = sqrt(N_current./N_prior);
pcer = 1/length(Z_prior);
%%
[~,~,~,q_min] = regularized_weights(mu,sigma,pcer);
%q_min is 0.2447

%% Data set 2:
prior = 'Lipids-TG-ONE.Teslovich.Nature.2010';
current = 'PGC-SCZ-2012.Nat.Genet.2012';
[~, ~, N_prior, Z_prior,~, N_current] = load_data(prior, current);
%%
mu = sqrt(N_current./N_prior).*Z_prior;
sigma = sqrt(N_current./N_prior);
pcer = 1/length(Z_prior);
%%
[~,~,~,q_min] = regularized_weights(mu,sigma,pcer);
%q_min is 0.2836

%% Data set 3:
prior = 'CKD-eGFRcrea.Kottgen.Nat.Genet.2010';
current = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';
[~, ~, N_prior, Z_prior,~, N_current] = load_data(prior, current);
%%
mu = sqrt(N_current./N_prior).*Z_prior;
sigma = sqrt(N_current./N_prior);
pcer = 1/length(Z_prior);
%%
[~,~,~,q_min] = regularized_weights(mu,sigma,pcer);
%q_min is 0.2337