function [w,num_loop_breaks, num_centering_steps] = monotone_weights_sub(mu,pcer,min_weight,max_weight,version,l,error_check)
% Compute optimal monotone weights for multiple testing.
% Use subsampling to compute things on a monotone grid.

% Inputs:
% mu - a negative vector of length J, the estimated means of test statistics
% pcer - per-comparison error rate at which the number of discoveries should be maximized;
%        pfer = J* pcer is the number of expected false rejections under the null
% min_weight - (optional) lower bound on weights, default  = 0;
% max_weight - (optional) upper bound on weights, default  = Inf;
% version - (optional) which version of subsampling should it use; default
% is 2
% error_check - (optional) display detailed error checking/diagnostics,
% default = 0


% Outputs:
% w - the optimal monotone weights that increase with the magnitude of mu

if ~exist('error_check','var')
    error_check = 0;
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

if ~exist('version','var')
    version = 2;
end

if ~exist('l','var')
    %maxi  = 5*1e3; %1e4
    maxi  = 1e4; %1e4
    l = min(maxi,length(mu));
end



%sort
[mu,ind] = sort(mu,'descend');
L=length(mu);

%subsample
switch version
    %version 1: keep a fixed number of means
    case 1
        ind_sub =floor(linspace(1,L,l));
        mu_s = mu(ind_sub);
        mu_s = sort(unique(mu_s),'descend');
        
        %version 2: keep all means above a fixed gap
    case 2
        %original
        ind_sub =floor(linspace(1,L,l));
        mu_s0 = mu(ind_sub);
        mu_s0 = sort(unique(mu_s0),'descend');
        L0=length(mu_s0);
        min_gap = 10^(-6);
        diff = mu_s0(1:L0-1)-mu_s0(2:L0);
        mu_s=mu_s0(L0);
        for i=1:L0-1
            if diff(L0-i)>min_gap
                mu_s = [mu_s0(L0-i); mu_s];
            end
        end
end
J = length(mu_s);

%strictly feasible starting point
%nearly equispaced weight spanning nearly the entire interval
w_sub = ((1:J)/(J+1))';
ep = max(0.9,(1-10/J))* min ( (max_weight-1)/(max(w_sub) - mean(w_sub)), (1 - min_weight)/(mean(w_sub)-min(w_sub)));
w_sub = 1 + ep*(w_sub-mean(w_sub));

if (any(w_sub>max_weight)||any(w_sub<min_weight))
    fprintf('Weights not in right interval\n');
end


t = 1e3;%in cvx book it states to choose this based on duality gap for some feasible var
mu_0 = 10; %based on cvx book '10-100 good'
epsi = 1e-5;
m = J+1; %number of constraints



%Barrier method iteration
num_centering_steps = 0;
num_loop_breaks = 0;
while m/t > epsi
    %Centering step
    [w_new,loop_break] = mon_weight_centering(mu_s,pcer,max_weight,min_weight,w_sub,t,epsi,error_check);
    num_loop_breaks = num_loop_breaks + loop_break;
    %Update
    w_sub = w_new;
    %increase t
    t = mu_0*t;
    num_centering_steps = num_centering_steps + 1;  
end

if (any(w_sub>max_weight)||any(w_sub<min_weight))
    fprintf('Weights not in right interval\n');
end
if any (w_new(2:J)-w_new(1:J-1) <0)
    fprintf('Weights not monotone\n');
end

%figure, plot(mu_sub,w_sub)

%extrapolate weights
%vq = interp1(x,v,xq) returns interpolated values of a 1-D function at specific query points using linear interpolation.
%Vector x contains the sample points, and v contains the corresponding values, v(x).
%Vector xq contains the coordinates of the query points.
w = interp1(mu_s,w_sub,mu,'spline');
%figure, plot(mu,w)

%sort weights according to mu
w_unsort = w;
for i=1:L
    w(ind(i))=w_unsort(i);
end

%ensure nonnegativity
%in some cases weights <0
w = max(w,0);
w = length(w)*w/sum(w); %renormalize