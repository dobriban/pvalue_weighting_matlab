function w = exp_weights(mu,beta,q)

J = length(mu);
s = zeros(J,1);

[mu_s,sort_index] = sort(mu);
u = exp(beta*mu_s);
c = mean(u);
S  = sum(s);
w = u/c*(J-S/q)/(J-S);
ind = find(w >1/q);

surplus = sum(w(ind))-length(w(ind))/q;
w(ind) = 1/q;
k = length(w(ind));

while (surplus>0)
increment = min(1/q-w(J-k),surplus);
w(J-k) = w(J-k) + increment;
surplus = surplus - increment;
k = k+1;
end

v = zeros(J,1);
for i=1:J
v(i) = w(find(sort_index==i));    
end
w = v;
%plot(mu,w,'*');

% eps = 1e-3;
% if abs(mean(w)-1)>eps
%     fprintf('Error: weights do not average 1');
% else
%     fprintf('OK\n');
% end