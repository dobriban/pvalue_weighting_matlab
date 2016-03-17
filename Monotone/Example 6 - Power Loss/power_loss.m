% evaluate power loss of monotone weights compared to optimal Spjotvoll
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 6 - Power Loss');
addpath '../../Code' '../../Code/Helper Code'

%%
savefigs =0; a = {'-','--','-.',':'};
J=1e4;
l_arr = [1e-3; 5*1e-3; 1e-2; 5*1e-2; 1e-1; 5*1e-1;0.9];
u_arr = [2; 10; Inf];
L = length(l_arr);
U = length(u_arr);
power_arr = zeros(L,U);
rng(0);  mu = -abs(randn(J,1));
pcer = 0.05/J;
%compute power if the truth is really fixed
power = @(w) sum(normcdf(norminv(pcer*w)-mu))/J;

[w_spjot,~] = spjotvoll_weights(mu,pcer);
p_spjot = power(w_spjot);
p_unw= power(ones(J,1));

for i=1:L
    l = l_arr(i);
    for j=1:U
        u = u_arr(j);
        tic
        w = monotone_weights_sub(mu,pcer,l,u);
        toc
        %find weights
        power_arr(i,j) = power(w);
    end
end
%%
%power_arr  = J*power_arr;
rng(2);
figure, hold on
h1 = plot(-log10(l_arr), power_arr(:,1),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(-log10(l_arr), power_arr(:,2),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(-log10(l_arr), power_arr(:,3),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(-log10(l_arr), ones(L,1)*p_unw,'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
h5 = plot(-log10(l_arr), ones(L,1)*p_spjot,'linewidth',4,'color',rand(1,3));
set(h5,'LineStyle',a{1});
n1 = num2str(u_arr);
legend([h1 h2 h3 h4 h5],{n1(1,:),n1(2,:),n1(3,:),'Unw','Spj'},'location','Best')
xlabel('-log_{10} l')
ylabel('Power');
xlim([min(-log10(l_arr)),max(-log10(l_arr))]);
set(gca,'fontsize',20)
if savefigs ==1
    filename = sprintf( './Img/power_loss_mon_weights_J_=_%d.png',J);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end