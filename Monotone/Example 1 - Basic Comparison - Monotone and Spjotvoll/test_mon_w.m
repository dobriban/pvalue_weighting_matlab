%test monotone weights
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 1 - Basic Comparison - Monotone and Spjotvoll');
addpath '../../Code' '../../Code/Helper Code'

%% basic computation of Spjotvoll & monotone
rng(0);
p = 1e3;
pcer = 5e-3;
mu = -abs(randn(p,1));
[w,c] =spjotvoll_weights(mu,pcer);
tic
[w_m] = monotone_weights_sub(mu,pcer);
toc
%% plot computed weights: linear and log-scale
a = {'-','--','-.',':'};
rng(2);
figure,
subplot(1,2,1),  hold on
[~,ind] = sort(mu);
h1=plot(mu(ind),w(ind),'*','linewidth',4,'color',rand(1,3));
%set(h1,'LineStyle','*');
h2=plot(mu(ind),w_m(ind),'-','linewidth',4,'color',rand(1,3));
%set(h2,'LineStyle','-');
legend('Spjotvoll','Monotone','Location','best');
xlabel('Prior Mean');
ylabel('Weight');
set(gca,'fontsize',20)
rng(2)
subplot(1,2,2),  hold on
h1=plot(mu(ind),log10(w(ind)),'*','linewidth',4,'color',rand(1,3));
%set(h1,'LineStyle',a{4});
h2=plot(mu(ind),log10(w_m(ind)),'-','linewidth',4,'color',rand(1,3));
%set(h2,'LineStyle',a{2});
legend('Spjotvoll','Monotone','Location','Southwest');
xlabel('Prior Mean');
ylabel('log10(Weight)');
set(gca,'fontsize',20)
%%
filename = sprintf( 'Monotone_and_Spjotvoll_Weights_Linear_Log.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
%% numerical comparison of Spjotvoll & monotone
%in cases where Spjotvoll is monotone, so they should agree within
%numerical precision
rng(0);
p = 1e3;
num_exp = 1e3;
beta = 1/2;
min_weight = 0;
err = zeros(num_exp,1);
corre = zeros(num_exp,1);
tolerance = 1e-2;
mu_matrix = zeros(num_exp,p);
weight_matrix_s = zeros(num_exp,p);
weight_matrix_m = zeros(num_exp,p);
large_error = zeros(num_exp,1);
q_vector = zeros(num_exp,1);

for i=1:num_exp
    mu = -sort(abs(randn(p,1)));
    mu_matrix(i,:) = mu;
    pcer = 5e-7*rand;
    [w,c] =spjotvoll_weights(mu,pcer);
    %ensure Spjotvoll monotone
    while any (w(2:p)-w(1:p-1) <0)
        pcer = beta*pcer;
        [w,c] =spjotvoll_weights(mu,pcer);
    end
    q_vector(i) = pcer;
    weight_matrix_s(i,:) = w;
    
    [w_m] = monotone_weights_sub(mu,pcer,min_weight);
    weight_matrix_m(i,:) = w_m;
    err(i) = sum(abs(w_m-w));
    corre(i) = corr(w_m,w);
end

%% Worst example
index = find(err==max(err));
mu  = mu_matrix(index,:);
w = weight_matrix_s(index,:);
w_m = weight_matrix_m(index,:);
a = {'-','--','-.',':'};
rng(2);
figure,  hold on
[~,ind] = sort(mu);
plot(mu(ind),w(ind),'*','linewidth',2,'color', rand(1,3))
plot(mu(ind),w_m(ind),'-','linewidth',2,'color', rand(1,3))
legend('Spjotvoll','Monotone','Location','Best');
xlabel('Prior Mean');
ylabel('Weight');
set(gca,'fontsize',20)
filename = sprintf( 'Numerical DisAgreement of Monotone and Spjotvoll.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Histograms of MAD & (1-corr)
figure,
subplot(1,2,1),  hold on
hist(log10(err/p),num_exp)
xlabel('log10(MAD)');
ylabel('Frequency');
set(gca,'fontsize',20)
subplot(1,2,2),  hold on
hist(log10(1-corre),num_exp)
xlabel('log10(1 - rho)');
ylabel('Frequency');
set(gca,'fontsize',20)
filename = sprintf( 'Numerical Agreement of Monotone and Spjotvoll.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Plot when both monotone
figure,  hold on
[~,ind] = sort(mu);
rng(2)
plot(mu(ind),w(ind),'*','linewidth',2,'color', rand(1,3))
plot(mu(ind),w_m(ind),'-','linewidth',2,'color', rand(1,3))
legend('Spjotvoll','Monotone','Location','best');
xlabel('Prior Mean');
ylabel('Weight');
set(gca,'fontsize',20)
filename = sprintf( 'Monotone_and_Spjotvoll_Both_Mon.png');
saveas(gcf, filename,'png');