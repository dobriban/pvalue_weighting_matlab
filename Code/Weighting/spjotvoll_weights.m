function [w,c,err] = spjotvoll_weights(mu,pcer)
% Compute optimal weights for multiple testing.
% Given estimated means of test statistics, the weighting scheme optimizes
% the expected number of discoveries at some specific level q
%
% Inputs:
% mu - a non-positive vector of length J, the estimated means of test statistics
% pcer - per-comparison error rate at which the number of discoveries should be maximized;
%        pfer = J* pcer is the number of expected false rejections under the null

% Outputs:
% w - the optimal weights. A non-negative vector of length J', where J' is
% the number of strictly negative means mu(i).
% c - the normalizing constant produced by solving the optimization problem.
%
% Author: Edgar Dobriban

err = 0;
J = length(mu);
pfer = J*pcer;
gamma = sum(normcdf(mu(mu<0)./2));
n0 = sum(mu==0);
neg_ind = (mu<0);
w = zeros(J,1);

%case 1: the mass at negative means is large enough
%then all weights are defined by the formula Phi(mu/2 + c/mu)
if (gamma >= pfer)
    [w_neg,c,err] = spjotvoll_weights_finder(mu(neg_ind),pcer,J);
    w(neg_ind) = w_neg;
else
    %case 2: the mass at negative means is small enough
    %then all weights are defined by the formula Phi(mu/2 + c/mu)
    if (gamma+n0 <= pfer)
        [w_neg,c,err] = spjotvoll_weights_finder(mu(neg_ind),J - n0/pcer);
        w(neg_ind) = w_neg;
        w(~neg_ind) = 1/pcer;
    else
        %case 3: the mass at negative means is intermediate
        w(neg_ind) = normcdf(mu(mu<0)./2)/pcer;
        w(~neg_ind) = (pfer-gamma)/(pcer*n0);
        c = NaN;
    end
end