function [delta_w_nt,lambda] = mon_weight_newton_step(mu,pcer,max_weight,min_weight,w,t,error_check)


J = length(mu); 
B = zeros(J,3);
q = pcer;
%set up columns of Hessian
%columns from log-penalty
l = 1./((w(2:J)-w(1:(J-1))).^2); %easily get numbers of order theta(J^2)
B((1:J-1),1) = -l;
B((1:J-1),2) = l;
B((2:J),2) = B((2:J),2) + l;
B((2:J),3) = -l;
B(1,2) = B(1,2) + 1./((w(1)-min_weight).^2);
B(J,2) = B(J,2) + 1./((w(J)-max_weight).^2);
%diagonal from main objective
B(:,2) = B(:,2) - sqrt(2*pi)*t.*mu.*q.^2.*exp(norminv(q.*w).^2./2 + mu.*norminv(q.*w)-mu.^2/2);
%form sparse tridiagonal mx (hessian of main objective)
T = spdiags(B,-1:1,J,J);

%solve Newton system & solve the two tridiagonal systems
grad_f = - t*q*exp(mu.*norminv(q.*w)-mu.^2/2) -1./[w(1)-min_weight; w(2:J)-w(1:(J-1))] + ...
    1./[w(2:J)-w(1:(J-1)); max_weight - w(J)];

d = ones(J,1);
%time sparse tridiagonal solve
tic 
s = T\d;
if (error_check==1)
toc
end
r = - 1/(d'*s);
u = -T\grad_f;
nu = r*(s'*grad_f);
delta_w_nt = u  - nu*s;

if (error_check==1)
fprintf('Newton step violates constraint by: %e\n',sum(delta_w_nt));
end

lambda = sqrt(delta_w_nt'*T*delta_w_nt);
