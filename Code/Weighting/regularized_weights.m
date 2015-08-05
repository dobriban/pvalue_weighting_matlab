function [w,lambda,true_constraint,q_min] = regularized_weights(mu,sigma,pcer,method)
% Compute the optimal regularized weights for p-value weighting.
% Given estimated means and standard errors of test statistics,
% the weighting scheme optimizes the expected number of discoveries
% at some specific level q
%
% Inputs:
% The inputs specify the distribution of the means as Normal(mu, sigma^2)
% where:
%
% mu - a negative vector of length J, the estimated means of test statistics
% sigma - a positive vector of length J, the estimated standard errors of test statistics
% pcer - per-comparison error rate at which the number of discoveries should be maximized;
%        pfer = J* pcer is the number of expected false rejections under the null
% method - (optional) 'Newton' (default) or 'Brent';  
%
% Outputs:
% w - the optimal weights. A non-negative vector of length J', where J' is
% the number of strictly negative means mu(i).
% lambda - the dual certificate
%   normalizing constant produced by solving the optimization problem.
% true_constraint - a variable indicating the true value of the constraint
% that was solved. This will equal J*pcer if the original problem was
% solved, or else it can equal a value close to it.
%
% Author: Edgar Dobriban



if (any(mu>=0))
  error('reg_w:negative_means', 'Negative means required')
end

if (any(sigma<=0))
  error('reg_w:positive_vars', 'Positive variances required')
end

sigma_threshold = 1e-2;
if (any(sigma<=sigma_threshold)) 
  sigma = sigma + sigma_threshold;
  fprintf('Warning: Some values of the standard error are small (below %f). Adding a small constant to them.\n',sigma_threshold);
end

if (pcer <=0) || (pcer >=1)
  error('reg_w:invalid_pcer', 'Per-comparison error rate must be in (0,1)')
end

J = length(mu);

if (length(sigma)==1)
  sigma = sigma*ones(J,1);
end

var_plus = (sigma.^2+1);

if ~exist('method','var')
  method = 'Newton';
end

switch method
  case 'Newton'
    
    var = sigma.^2;
    mu2 = mu.^2;
    alpha = mu./var;
    beta = (mu2 + var.*(mu2+var_plus.*log(var_plus)))./(var.^2);
    gamma = 2*var_plus./var;
    
    %the function c_1 in my paper equals
    %c_1(eta,sigma,lambda) = -alpha - sqrt(beta + gamma.*c)
    % where c = log(lambda)
    %so checking the condition for dual feasibility corresponds to verifying
    % f(0) >= 0
    %for f defined below
    f = @(c) 1/J*sum(normcdf(-alpha - sqrt(beta + gamma.*c)))-pcer;
    df = @(c) 1/J*sum(-normpdf(-alpha - sqrt(beta + gamma.*c)).*gamma.*(beta + gamma.*c).^(-1/2));
    
    if (f(0)<0)
      error('reg_w:dual_var', 'Dual certifiability condition failed')
    end
    %set the value of minimal q
    q_min = f(0)+pcer;
    
    x0 = 0;
    epsi = 1e-3;
    tol = epsi*pcer/J;
    nmax = 100;
    
    % code from Matlab file-exchange function 'newton'
    % Author:	Tashi Ravach
    x(1) = x0 - (f(x0)/df(x0));
    ex(1) = abs(x(1)-x0);
    k = 2;
%     while (abs(ex(k-1)) >= tol) && (k <= nmax)
%       x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
%       ex(k) = abs(x(k)-x(k-1));
%       k = k+1;
%     end
    
    x(1) = x0 - (f(x0)/df(x0));
    k = 2;
    while (abs(f(x(k-1))) >= tol) && (k <= nmax)
      x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
      k = k+1;
    end
    %plot function values with iteration
    %     f_x = zeros(k-1,1);
    %     for i=1:k-1
    %     f_x(i) = f(x(i));
    %     end
    %     figure, plot((1:k-1),f_x,'*');
    %
    lambda = exp(x(k-1));
    
    w = normcdf(-alpha - sqrt(beta + gamma.*x(k-1)))/pcer;
    true_constraint = J*pcer;
  case 'Brent'
    eta = mu; %this is just a naming convention
    gamma = sqrt(var_plus);
    q = pcer;
    l_prime = zeros(J,1);
    for i=1:J
      l_prime(i) = find_crossing(eta(i),gamma(i));
    end
    
    %find the interval where the right dual variable may lie
    [l_sort,ind] = sort(l_prime);
    %l_sort = l_prime(ind)
    %ind(k) is the index of the k-th smallest one
    
    %global dual function
    H = @(lambda) sum(g_fun(eta,gamma,lambda,l_prime));
    
    %if the sought-for interval is between the l_prime values
    %artificially augment l_sort with a 1
    l_sort  = [l_sort; 1];
    if H(l_sort(J+1))<=J*q
      %binary search: l - lower; u - upper;
      l=1; u=J+1;
      while (u-l>1)
        mid = floor((u+l)/2);
        if H(l_sort(mid))<J*q
          u=mid;
        else
          l=mid;
        end
      end
      %right interval is l,l+1
      a = H(l_sort(l));
      c = H(l_sort(l+1));
      
      %check a <= J;
      %this should never happen, because it corresponds to lambda > 1
      if a>J
        fprintf('Error: dual constraint violated at 1');
      end
      
      %right limit of function at l_sort(l)
      b = a - 1 + normcdf(c_1(eta(ind(l)),gamma(ind(l)),l_sort(l)));
      
      %should have a>b>c
      if (a<b)||(b<c)
        fprintf('Error: incorrect ordering of interval endpoints');
      end
      if (J*q>a)||(J*q<c)
        fprintf('Error: interval does not contain target');
      end
      
      % so we know c < b < a, and that c < Jq < a
      %two cases, depending on the ordering of Jq, b
      %(1) ; Jq is in (b,a): in this case can't certify the original problem
      if (J*q > b)
        lambda = l_sort(l);
        true_constraint = a; %instead of Jq
        true_q = a/J;
        
        %(2) : Jq is in (c,b): in this case can certify the original problem
      else
        %find the zero of the following function f
        f = @(lambda) H(lambda) - J*q;
        %will use standard matlab zero-finding; need starting interval
        %ensure starting points have opposite sign: f is decreasing
        epsi = 1e-7;
        l_1 = l_sort(l) + epsi;
        l_2 = l_sort(l+1);
        while(f(l_1) < 0)
          epsi = epsi/10;
          l_1 = l_sort(l) + epsi;
        end
        x0 = [l_1 l_2]; % initial interval
        
        lambda = fzero(f,x0);
        true_constraint = J*q;
        true_q = q;
      end
      
    else %if the right interval is bigger than 1
      f = @(lambda) H(lambda) - J*q;
      %ensure starting points have opposite sign
      epsi = 1e-7;
      m = 1;
      l_1 = m + epsi;
      l_2 = 1;
      while(f(l_1) < 0)
        epsi = epsi/10;
        l_1 =m + epsi;
      end
      big_const = 1;
      while(f(l_2) > 0)
        big_const = big_const*2;
        l_2 =m + big_const;
      end
      x0 = [l_1 l_2]; % initial interval
      
      lambda = fzero(f,x0);
      true_constraint = J*q;
      true_q = q;
    end
    
    %calculate weights
    %note: dividinf by true_q instead of q
    w =  g_fun(eta,gamma,lambda,l_prime)/true_q;
end

%error check:
%epsi = tol*J;
  
if abs(sum(w)-J)>epsi
  fprintf('Warning: Weights do not average to 1. Mean weight - 1 = %E\n',mean(w)-1);
  %error('reg_w:invalid_weights_sum', 'Weights do not average 1')
end
if any(w<0)
  error('reg_w:invalid_weights_neg', 'Some weights are negative')
end
