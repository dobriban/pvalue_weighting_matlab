%compare with R
J = 200;
mu = -abs(randn(J,1));
sigma = ones(J,1);
gamma = sqrt(1+sigma.^2);
l_prime = zeros(J,1);
for i=1:J
    l_prime(i) = find_crossing(mu(i),gamma(i));
end
plot(mu,l_prime,'*')

%%
lambda = 0.58;
g = g_fun(mu,gamma,lambda,l_prime);
plot(mu,g,'*')
%%
alpha = 100;
w_1 = bayes_weights(mu, sigma, alpha/J);
plot(mu,w_1,'*')

%%
eta = -1;
gamma = 2;
lambda_grid = (1:J)'/J;
for i=1:J
g(i) = g_fun(eta,gamma,lambda_grid(i));
end
plot(lambda_grid,g)