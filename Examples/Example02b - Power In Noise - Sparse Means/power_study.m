%% Study the power of various weighting schemes
% depends on weighting_power_study
codedir =  ['C:/Dropbox/Weighted New/pvalue_weighting_matlab/Examples/Example02b - Power In Noise - Sparse Means'];
cd(codedir);
%
addpath '../../Code'
addpath '../../Code/Helper Code'
%%
J = 1e3; %number of tests
%sparsity fraction
sparsity = linspace(0.01,0.1,10);
s = length(sparsity);
pfer = 10; %expected number of false rejections under 'null'
m = 1e-3; %size of small mean
M = 2; %size of large mean
sigma = 1;

pd_unif = zeros(s,1);
pd_spjot = zeros(s,1);
pd_reg = zeros(s,1);
pa_unif = zeros(s,1);
pa_spjot = zeros(s,1);
pa_reg = zeros(s,1);

w_s = zeros(s,2);
w_r = zeros(s,2);

for i=1:s
  pi_0 = 1-sparsity(i);
  [pd_unif(i), pd_spjot(i), pd_reg(i), pa_unif(i), pa_spjot(i), pa_reg(i), w_spjot,w_reg] = ...
    power_study_helper(J,pfer,pi_0,m,M,sigma);
  %save weights
  w_s(i,:) = [w_spjot(1), w_spjot(J)];
  w_r(i,:) = [w_reg(1), w_reg(J)];
end
%%
a = {'-','--',':','-.'};
figure
subplot(1,2,1)
set(gca,'fontsize',8)
col = colormap(gray(4));

h = plot(sparsity, pd_unif); hold on
set(h,'LineStyle',a{1});
set(h,'LineWidth',2);
%set(h,'color',col(1,:));
h = plot(sparsity, pd_spjot,'g'); hold on
set(h,'LineStyle',a{2});
set(h,'LineWidth',2);
%set(h,'color',col(2,:));
h = plot(sparsity, pd_reg,'r'); hold on
set(h,'LineStyle',a{3});
set(h,'LineWidth',2);
%set(h,'color',col(3,:));
%plot(sparsity, pd_unif,sparsity, pd_spjot,'g-.*',sparsity, pd_reg,'r-o')
title('Deterministic');
xlabel('\pi_1');
hleg1 = legend('Unw','Spjot', 'Info');
set(hleg1,'Location','North')
set(hleg1,'FontSize',7);
mi = [min(pd_unif),min(pd_spjot),min(pd_reg)];
ma = [max(pd_unif),max(pd_spjot),max(pd_reg)];
ylim([min(mi),max(ma)*1.4])
set(gca,'fontsize',8)

subplot(1,2,2)
h = plot(sparsity, pa_unif); hold on

set(h,'LineStyle',a{1});
set(h,'LineWidth',2);
h = plot(sparsity, pa_spjot,'g'); hold on
set(h,'LineStyle',a{2});
set(h,'LineWidth',2);
h = plot(sparsity, pa_reg,'r'); hold on
set(h,'LineStyle',a{3});
set(h,'LineWidth',2);
%plot(sparsity, pa_unif,sparsity, pa_spjot,'g-.*',sparsity, pa_reg,'r-o')
title('Average');
xlabel('\pi_1');
set(gca,'fontsize',8)
%%
saveTightFigure(gcf,'PowerComparison_SparseMeans.pdf')
%% plot the weights themselves
a = {'-','--',':','-.'};
figure,
subplot(1,2,1)
h = plot(sparsity, w_s(:,2)); hold on
set(h,'LineStyle',a{1});
set(h,'LineWidth',2);
h =plot(sparsity, w_s(:,1));
set(h,'LineWidth',2);
set(h,'LineStyle',a{2});
title('Spjotvoll');
xlabel('\pi_1');
ylim([min(min(w_s)),max(max(w_s))*1.1])
set(gca,'fontsize',8)

subplot(1,2,2)

h = plot(sparsity, w_r(:,2)); hold on
set(h,'LineStyle',a{1});
set(h,'LineWidth',2);
h =plot(sparsity, w_r(:,1));
set(h,'LineWidth',2);
set(h,'LineStyle',a{2});
title('Informed');
xlabel('\pi_1');
hleg1 = legend('Large', 'Small');
set(hleg1,'Location','NorthWest')
set(hleg1,'FontSize',7);
ylim([min(min(w_r)),max(max(w_r))*1.4])
set(gca,'fontsize',8)


saveTightFigure(gcf,'PowerComparison_SparseMeans_Weights.pdf')