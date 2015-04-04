function [pd_unif, pd_spjot, pd_regularized, pa_unif, pa_spjot, pa_regularized] ...
  = power_study_helper(J,pfer,c,noise)
%helper function in the power study of p-value weighting

% Inputs:
% J - number of means to be generated
% pfer - per-family error rate at which the number of discoveries should be maximized; 
%        pfer is the number of expected false rejections under the null
% c - sparsity size, a positive real >0  
%     means are generated proportional to betarnd(c,1)
% noise - noise level to be used in the experiment

%Outputs: Powers of three procedures: Uniform, Spjotvoll and Regularized.
%pd_ - power of the procedures if the means are truly equal to the value
%that the weighting scheme uses
%pa_ - power of the procedures if the means are random

epsi = 1e-3;
grid = linspace (epsi,1-epsi,J)';

mu = - 4*betainv(grid,c,1);
sigma = noise.*ones(J,1); 
alpha = pfer/J;

%find weights
[w_spjot,~] = spjotvoll_weights(mu,alpha);  
[w_reg,~] = bayes_weights(mu,sigma,alpha);

%compute power if the truth is really fixed
power = @(w) sum(normcdf(norminv(alpha*w)-mu))/J;
power0 = @(w) sum(normcdf(norminv(alpha*w)-mu(mu<0)))/J;
pd_unif = power(ones(J,1)); %power_deterministic
pd_spjot = power0(w_spjot);
pd_regularized = power0(w_reg);

%compute power if the truth is really random
power = @(w) sum(normcdf((norminv(alpha*w)-mu)./sqrt(sigma.^2+1)))/J;
power0 = @(w) sum(normcdf((norminv(alpha*w)-mu(mu<0))./sqrt(sigma(mu<0).^2+1)))/J;
pa_unif = power(ones(J,1)); %power_average
pa_spjot = power0(w_spjot);
pa_regularized = power0(w_reg);
