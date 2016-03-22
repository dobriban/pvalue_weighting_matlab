function w = exp_weights_unbounded(mu,beta)
%remove bound from exp_weights

u = exp(-beta*mu);
w = u/mean(u);
