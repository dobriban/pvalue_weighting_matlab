function [lambda] = find_crossing(eta,gamma,plot)
%find the crossing point of the two competing maxima
%Inputs:
%eta - mean
%gamma > 1 - marginal variance: gamma^2 = 1 + sigma^2
%
%lambda - crossing point of the two maxima
% work with x = log(lambda)

%correspondence between variables 
%lambda:  - l
%x          : f log_l 


f = @(x) 1  - exp(x) - (normcdf((c_1(eta,gamma,x,1) - eta)/gamma) - ...
  exp(x)*normcdf((c_1(eta,gamma,x,1))));

%f is a decreasing function that has a zero in [log_l, 0]
%this is log(l), where l is the lower bound on the crossing in lambda-space
log_l =  - log(gamma)  - eta^2/(2*(gamma^2-1));
%in principle the starting interval is [log_l,0]

%initialize to ensure that starting point obeys:
%square(x) >0
%f(x) > 0
% square = @(x) eta_scalar^2 + (gamma_scalar^2-1)*(eta_scalar^2 + 2*gamma_scalar^2*(log(gamma_scalar) + x));
% epsi = 1e-3;
% 
% 
% while (square(log_l) <0)
%   while(f(log_l + epsi) < 0)
%     epsi = epsi/10;
%   end
%   log_l = log_l + epsi;
% end

if (f(log_l) > 0)
  x0 = [log_l 0]; % initial interval
  x = fzero(f,x0);
  
  lambda = exp(x);
else %within numerical precision, lambda equals the starting point
  lambda = exp(log_l);
end
%% optional plotting feature;
if exist('plot','var')
  
  nil = @(x) 0;
  figure
  hold on
  fplot(f,x0)
  fplot(nil,x0)
  plot(x,0,'ro')
  
end