%% Study the power of various weighting schemes
% depends on weighting_power_study
codedir =  ['C:/Git/pvalue_weighting_matlab/Examples/Example10- Compare Weighting Methods'];
cd(codedir);
addpath '../../Code'
addpath '../../Code/Helper Code'
%% low noise absolute normal sigma

J = 1e3; %number of tests
rng(0);
sigma = abs(randn(J,1));
mu = randn(J,1);
pfer = 10;
pcer = pfer/J;
q = pcer;

%power0 = @(w) sum(normcdf((norminv(pcer*w(mu<0))-mu(mu<0))./sqrt(sigma(mu<0).^2+1)))/J;
power0 = @(w) sum(normcdf((norminv(min(pcer*w,1))-mu)./sqrt(sigma.^2+1)))/J;

global_array = linspace(0.01,4,100);
dispersion = global_array;
for i=1:length(dispersion);
    w = bayes_weights(mu,dispersion(i)*sigma,pcer);
    power_reg(i) = power0(w);
end

beta = global_array;
for i=1:length(beta);
    w = exp_weights(abs(mu), beta(i), pcer);
    w = w/mean(w);
    power_exp(i) = power0(w);
end

threshold = global_array;
[mu_s,sort_index] = sort(mu);
k = floor(J*q);
topw = sort_index(1:k);

for i=1:length(threshold);
    ind = find(mu<-threshold(i));
    w = zeros(J,1);
    w(ind) = 1;
    if mean(w)>0
        w = w/mean(w);
        if (max(w) > 1/q)
            w = zeros(J,1);
            w(topw) = 1/q;
        end
        power_top(i) = power0(w);
    else
        w = zeros(J,1);
        w(topw) = 1/q;
        power_top(i) = power0(w);  
    end
end

power_unw = power0(ones(J,1))*ones(length(global_array),1);
%%
a = {'-','--',':','-.'};
col = colormap(gray(5));
figure, hold on
h = plot(global_array, power_unw,'linewidth',1); set(h,'LineStyle',a{1}); set(h,'color',col(1,:));
h = plot(global_array, power_reg,'linewidth',1); set(h,'LineStyle',a{2});set(h,'color',col(2,:));
h = plot(global_array, power_exp,'linewidth',1); set(h,'LineStyle',a{3});set(h,'color',col(3,:));
h = plot(global_array, power_top,'linewidth',1); set(h,'LineStyle',a{4});set(h,'color',col(4,:));
ylim([0,0.12])
%h = legend('Unweighted','Informed','Exponential','Filtering','location','Best');
ylabel('Power');
%set(h,'FontSize',8);

%%
saveTightFigure(gcf,'compare_weighting_methods_grayscale.pdf')
%% plot
filename = sprintf( './compare_weighting_methods.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);