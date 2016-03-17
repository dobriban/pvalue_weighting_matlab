function [w,num_loop_breaks, num_centering_steps, w_array] = monotone_weights(mu,pcer,min_weight,max_weight,error_check, all_iter)
% Compute optimal monotone weights for multiple testing.
% Inputs:
% mu - a negative vector of length J, the estimated means of test statistics
% pcer - per-comparison error rate at which the number of discoveries should be maximized;
%        pfer = J* pcer is the number of expected false rejections under the null
% min_weight - (optional) lower bound on weights, default  = 0;
% max_weight - (optional) upper bound on weights, default  = Inf;
% error_check - (optional) display detailed error checking/diagnostics,
% default = 0

% Outputs:
% w - the optimal monotone weights that increase with the magnitude of mu

if ~exist('error_check','var')
    error_check = 0;
end

if ~exist('all_iter','var')
    all_iter = 0;
end

if (any(mu>=0))
    error('reg_w:negative_means', 'Negative means required')
end

if ~exist('min_weight','var')
    min_weight = 0;
end

if ~exist('max_weight','var')
    max_weight = Inf;
end
max_weight = min(max_weight,1/pcer);

[mu,ind] = sort(mu,'descend');

%strictly feasible starting point
J=length(mu);
%nearly equispaced weight spanning nearly the entire interval
w = ((1:J)/(J+1))';
ep = max(0.9,(1-10/J))* min ( (max_weight-1)/(max(w) - mean(w)), (1 - min_weight)/(mean(w)-min(w)));
w = 1 + ep*(w-mean(w));

if (any(w>max_weight)||any(w<min_weight))
    fprintf('Weights not in right interval\n');
end


t = 1e3;%in cvx book it states to choose this based on duality gap for some feasible var
mu_0 = 10; %based on cvx book '10-100 good'
epsi = 1e-5;

m = J+1; %number of constraints


if (all_iter==1)
    w_array  = w;
end

%Barrier method iteration
num_centering_steps = 0;
num_loop_breaks = 0;
while m/t > epsi
    %Centering step
    [w_new,loop_break] = mon_weight_centering(mu,pcer,max_weight,min_weight,w,t,epsi,error_check);
    num_loop_breaks = num_loop_breaks + loop_break;
    %Update
    w = w_new;
    %increase t
    t = mu_0*t;
    num_centering_steps = num_centering_steps + 1;
    
    if (all_iter==1)
        w_array  = [w_array w];
    end
end

if (any(w>max_weight)||any(w<min_weight))
    fprintf('Weights not in right interval\n');
end
if any (w_new(2:J)-w_new(1:J-1) <0)
    fprintf('Weights not monotone\n');
end


%sort weights according to mu
w_unsort = w; 
for i=1:J
    w(ind(i))=w_unsort(i);
end

if (all_iter==1)
    w_array_unsort = w_array;
    for i=1:J
        w_array(ind(i),:)=w_array_unsort(i,:);
    end
end