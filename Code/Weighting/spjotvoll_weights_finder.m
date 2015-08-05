function [w,lambda,err] = spjotvoll_weights_finder(mu,pcer,t)
%helper function in finding weights
%maximize sum P(P_i \le  pcer * w_i)
%subject to sum w_i = t

err = 0;
J = length(mu);

f = @(c) 1/J*sum(normcdf(c./mu+mu/2))-pcer*t/J;
df = @(c) 1/J*sum(normpdf(c./mu+mu/2)./mu);

x0 = 0;
epsi = 1e-3;
tol = epsi*pcer/J;
nmax = 100;

% solve f = 0 using Newton's method
x(1) = x0 - (f(x0)/df(x0));
k = 2;
while (abs(f(x(k-1))) >= tol) && (k <= nmax)
    x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
    k = k+1;
end
% plot function values with iteration
% f_x = zeros(k-1,1);
% for i=1:k-1
%     f_x(i) = f(x(i));
% end
% figure, plot((1:k-1),log(abs(f_x)),'*');


c = x(k-1);
w = normcdf(c./mu+mu/2)/pcer;

lambda = exp(x(k-1));

if abs(sum(w)-t)>epsi
  fprintf('Warning: Weights do not average to 1. Mean weight - 1 = %E\n',mean(w)-1);
  err = 1;
  %error('reg_w:invalid_weights_sum', 'Weights do not average 1')
end
if any(w<0)
  error('reg_w:invalid_weights_neg', 'Some weights are negative')

end
