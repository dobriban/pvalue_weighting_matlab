function [pd_unif, pd_spjot, pd_regularized, pa_unif, pa_spjot, pa_regularized, w_spjot,w_reg] ...
  = power_study_helper(J,pfer,pi_0,m,M,sigma)
%helper function in the power study of p-value weighting
%use the two-point mixture pi_0 * delta_m + (1-pi_0)* delta_M

% Inputs:
% J     - number of means to be generated
% pfer  - per-family error rate at which the number of discoveries should be maximized; 
%        pfer is the number of expected false rejections under the null
% pi_0  - proportion of small means m
% m     - size of the small means
% M     - size of the large means
% sigma - noise level to be used in the experiment

%Outputs: Powers of three procedures: Uniform, Spjotvoll and Regularized.
%pd_ - power of the procedures if the means are truly equal to the value
%that the weighting scheme uses
%pa_ - power of the procedures if the means are random


n = floor(pi_0*J); %number of small means
mu = - [m*ones(n,1); M*ones(J-n,1)];
sigmav = sigma*ones(J,1); 
alpha = pfer/J;

%find weights
[w_spjot,~] = spjotvoll_weights(mu,alpha);
[w_reg,~] = bayes_weights(mu,sigmav,alpha); 

%compute power if the truth is really fixed
power = @(w) sum(normcdf(norminv(alpha*w)-mu))/J;
power0 = @(w) sum(normcdf(norminv(alpha*w)-mu(mu<0)))/J;
pd_unif = power(ones(J,1)); %power_deterministic
pd_spjot = power0(w_spjot);
pd_regularized = power0(w_reg);

%compute power if the truth is really random
power = @(w) sum(normcdf((norminv(alpha*w)-mu)./sqrt(sigmav.^2+1)))/J;
power0 = @(w) sum(normcdf((norminv(alpha*w)-mu(mu<0))./sqrt(sigmav(mu<0).^2+1)))/J;
pa_unif = power(ones(J,1)); %power_average
pa_spjot = power0(w_spjot);
pa_regularized = power0(w_reg);

