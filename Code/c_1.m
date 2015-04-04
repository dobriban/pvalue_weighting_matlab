function c1 = c_1(eta,gamma,lambda,logscale)
%the c_1 function gives the location of the smaller maximum of the generic
%term in the dual
%Inputs
%eta, gamma - prior mean and standard error
%lambda - the dual variable
%logscale - (optional) binary variable indicating if log-scale should be used
%for lambda
%      - default -0 -  no

if ~exist('logscale','var')
    logscale = 0;
end

%not log-scale (original version)
if (logscale == 0)
    if ~(lambda>=0)
        fprintf('lambda = %f\n',lambda);
        error('c_1:neg_lambda', 'Lambda is not positive')
    end
    x = log(lambda);
    
    %log-scale
else
    %input lambda is on log-scale
    x = lambda;
end

square = eta^2 + (gamma^2-1)*(eta^2 + 2*gamma^2*(log(gamma) + x));

%numerical slack
epsi = 1e-7;
if ~(square>=-epsi)
    fprintf('square = %E\n',square);
    error('c_1:neg_square', 'Quadratic expression q is not positive')
end
c1  =  - (eta + sqrt(abs(square)))/(gamma^2 - 1);
