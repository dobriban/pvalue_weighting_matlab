function pow = power_opt(mu_vec,pi_1_vec,q)
%calculate power of optimal weighting procedure
%when population distribution if a mixture of a point mass at 0 and at mu
%with frequencies pi_0 and pi_1

l1 = size(mu_vec,1);
l2 = size(mu_vec,2);

pow = zeros(l1,l2);

for i=1:l1
  for j=1:l2
  mu = mu_vec(i,j);
  pi_1 = pi_1_vec(i,j);
  
  if (pi_1 * normcdf(mu/2) > q)
    pow(i,j) = pi_1.*normcdf(norminv( q./pi_1) - mu);
  else
    pow(i,j) = q + pi_1.*(normcdf(- mu./2)-normcdf(mu./2));
  end
  
  end
end