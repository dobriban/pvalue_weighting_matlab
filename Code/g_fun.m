function [g,l_prime] = g_fun(eta,gamma,lambda,l_prime)
%g is the dual constraint function: 1 if lambda <= l_prime, and a phi(c_1)
%otherwise
%
%Inputs
% all are supposed to be of the same length
% eta, gamma - prior mean and standard error
% l_prime -(optional) precomputed crossing point
%
%lambda - the dual variable

J  = length(eta);
%compute l_prime on the fly if needed
if ~exist('l_prime','var')
  l_prime = zeros(J,1);
  for i=1:J
    [l_prime(i)] = find_crossing(eta(i),gamma(i));
  end
end

g = zeros(J,1);
for i=1:J
    if lambda <=l_prime(i)
        g(i) = 1;
    else
        g(i) = normcdf(real(c_1(eta(i),gamma(i),lambda)));
    end
end
