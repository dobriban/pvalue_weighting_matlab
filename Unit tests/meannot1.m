%Author: Joanna Bogdan
%Date: 5/18/2016
%Modified by Edgar Dobriban

%Error about not summing to J for situations with small sparsity (big k,
%depending on q as well
%doesnt sum to J for example for q=0.5,k=5 (does for q=0.5,k=3)
%doesnt for q=500,k=2 (does for q=500,k=1)


J = 1e3;
% q=0.0005;
% k=5;% varies sparsity - the bigger the k the less alternatives
q=0.5;
k=5;% varies sparsity - the bigger the k the less alternatives
rng(2)
eta = k + randn(J,1); %normal with mean m, standard deviation s (may be changed for different sparsity)
sigma=abs(randn(J,1));
alpha=q/J;
%%
[w,q_star,q_threshold, lambda]  = bayes_weights(eta,sigma,alpha);
err = abs(mean(w)-1);

%% Notes
%In bayes_weights, the values of l_prime are all nearly 1;
%l_prime = find_crossing(mu,gamma)
%length(unique(l_prime)) = 234;
%add jitter

%then in bayes_weights, the normalization of q in the last branch of the
%search for lambda was not correct
%I have now changed it to 
%q_star = Jq_temp;
%true_q = q_star/J;
        