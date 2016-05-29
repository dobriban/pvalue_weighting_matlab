%Author: Joanna Bogdan
%Date: 5/18/2016

%Error about not summing to J for situations with small sparsity (big k,
%depending on q as well

%doesnt sum to J for example for q=0.5,k=5 (does for q=0.5,k=3)
%doesnt for q=500,k=2 (does for q=500,k=1)


J = 1e3; 
q=0.0005;
k=5;% varies sparsity - the bigger the k the less alternatives
s=1;
eta = k + s.*randn(J,1); %normal with mean m, standard deviation s (may be changed for different sparsity)
sigma=abs(randn(J,1));
mu = normrnd(eta,sigma);
X = normrnd(mu,1);
P_current = normcdf(X); %current p-values;
alpha=q/J;

  [w,q_star,q_threshold, lambda]  = bayes_weights(eta,sigma,alpha);
  sum(w)
  
  