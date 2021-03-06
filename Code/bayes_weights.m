function [w,q_star,q_threshold, lambda] = bayes_weights(mu,sigma,q)
% Compute the optimal Bayes p-value weights
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
% q -  corrected bonferroni level (e.g. 0.05/J)
%
% Outputs:
% w - the optimal weights
% lambda - the dual certificate, normalizing constant produced by solving the optimization problem.
% q_star - true value of q solved for
% q_threshold - maximal value of q for which problem can be solved exactly
%
% Author: Edgar Dobriban



% if (any(mu>=0))
%   error('reg_w:negative_means', 'Negative means required')
% end

if (any(sigma<=0))
    error('reg_w:positive_vars', 'Positive variances required')
end

sigma_threshold = 1e-2;
if (any(sigma<=sigma_threshold))
    sigma = sigma + sigma_threshold;
    fprintf('Warning: Some values of the standard error are small (below %f). Adding a small constant to them.\n',sigma_threshold);
end

if (q <=0) || (q >=1)
    error('reg_w:invalid_pcer', 'Per-comparison error rate must be in (0,1)')
end

J = length(mu);

if (length(sigma)==1)
    sigma = sigma*ones(J,1);
end

var_plus = (sigma.^2+1);
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
f = @(c) 1/J*sum(normcdf(-alpha - sqrt(beta + gamma.*c)))-q;
df = @(c) 1/J*sum(-normpdf(-alpha - sqrt(beta + gamma.*c)).*gamma.*(beta + gamma.*c).^(-1/2));
q_threshold = f(0)+q;

if (f(0)>=0)
    x0 = 0;
    epsi = 1e-3;
    tol = epsi*q/J;
    nmax = 100;
    
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
    
    w = normcdf(-alpha - sqrt(beta + gamma.*x(k-1)))/q;
    q_star = q;
    
else %brent
    epsi = 1e-3;
    gamma = sqrt(var_plus);
    l_prime = zeros(J,1);
    overestimate_ind=0;
    for i=1:J
        l_prime(i) = find_crossing(mu(i),gamma(i));
    end
    
    %add some small jitter to deal with equal l-values
    if (length(unique(l_prime))<J)
        jitter_sd = 1e-2;
        l_jitter = l_prime;
        while (length(unique(l_jitter))<J)
            l_jitter = l_prime + jitter_sd*randn(J,1);
            jitter_sd = 2*jitter_sd;
        end
        l_prime = l_jitter;
    end
    %find the interval where the right dual variable may lie
    [l_sort,ind] = sort(l_prime);
    %l_sort = l_prime(ind)
    %ind(k) is the index of the k-th smallest one
    
    %global dual function
    H = @(lambda) sum(g_fun(mu,gamma,lambda,l_prime));
    
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
        b = a - 1 + normcdf(c_1(mu(ind(l)),gamma(ind(l)),l_sort(l)));
        
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
            %find closest endpoint
            dist_up = a - J*q;
            dist_down = J*q-b;
            if dist_up<dist_down
                q_star = a; %instead of Jq
                true_q = a/J;
            else
                q_star = b; %instead of Jq
                true_q = b/J;
                overestimate = (a-b);
                overestimate_ind = ind(l);
            end
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
            q_star = J*q;
            true_q = q;
        end
        
    else %if the right interval is bigger than 1, i.e., [1, +\infty)
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
        big_const = 1.01;
        while(f(l_2) > 0)
            l_2 =l_2*big_const;
        end
        x0 = [l_1 l_2]; % initial interval
        
        lambda = fzero(f,x0);
        Jq_temp = H(lambda); %should equal J*q
        %J x q_star = Jq_temp
        q_star = Jq_temp;
        true_q = q_star/J;
        
        %Old:
        %         q_star = J*q;
        %         true_q = q;
        %ezplot(f,x0);
    end
    
    %calculate weights
    %note: dividing by true_q instead of q
    w_temp =  g_fun(mu,gamma,lambda,l_prime);
    if overestimate_ind>0
        w_temp(overestimate_ind) = w_temp(overestimate_ind) - overestimate;
    end
    w =  w_temp/true_q;
end

%error check:
if abs(sum(w)-J)>epsi
    fprintf('Warning: Weights do not average to 1. Mean weight - 1 = %E\n',mean(w)-1);
    %error('reg_w:invalid_weights_sum', 'Weights do not average 1')
end
if any(w<0)
    error('reg_w:invalid_weights_neg', 'Some weights are negative')
end
