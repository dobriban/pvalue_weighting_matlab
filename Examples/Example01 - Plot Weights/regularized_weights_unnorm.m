function [w] = regularized_weights_unnorm(mu,sigma,pcer)
%generate un-normalized regularized weights
%used for plotting purposes

var_plus = (sigma.^2+1);
var = sigma.^2;
mu2 = mu.^2;
alpha = mu./var;
beta = (mu2 + var.*(mu2+var_plus.*log(var_plus)))./(var.^2);
gamma = 2*var_plus./var;

c = 2;

w = normcdf(-alpha - sqrt(beta + gamma.*c))/pcer;