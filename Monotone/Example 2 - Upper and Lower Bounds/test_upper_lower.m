%test monotone weights: upper and lower bounds
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 2 - Upper and Lower Bounds');
addpath '../../Code' '../../Code/Helper Code'

%% test changing the minimum weight
rng(0);
pcer = 5e-2;
p = 1e3;
mu = -sort(abs(randn(p,1)),'ascend');

%% plot all in one figure
M = 4;%used to be 10
min_weights =((0:(M-1))/M)';
w_m = zeros(length(min_weights),p);
for i=1:length(min_weights)
    min_weight  = min_weights(i);
    [w_m(i,:)] = monotone_weights_sub(mu,pcer,min_weights(i));
end
a = {'-','--','-.',':'};
rng(2);
figure, hold on
h1 = plot(mu,w_m(1,:),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(mu,w_m(2,:),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(mu,w_m(3,:),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(mu,w_m(4,:),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
xlabel('Prior Mean');
ylabel('Weight');
xlim([min(mu),max(mu)]);
set(gca,'fontsize',20)
filename = sprintf( 'Lower_Bound_weights_.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

%% test changing the maximum weight
rng(0);
pcer = 5e-2;
p = 1e3;
mu = -sort(abs(randn(p,1)),'ascend');
l = 1e3;

max_weights = 1 + (1:1:M)'/M;
w_m = zeros(length(max_weights),p);
for i=1:length(max_weights)
    max_weight  = max_weights(i);
    min_weight = 0;
    [w_m(i,:)] = monotone_weights_sub(mu,pcer,min_weight,max_weight);
end

rng(2);
figure, hold on
h1 = plot(mu,w_m(1,:),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(mu,w_m(2,:),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(mu,w_m(3,:),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(mu,w_m(4,:),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
xlabel('Prior Mean');
ylabel('Weight');
xlim([min(mu),max(mu)]);
set(gca,'fontsize',20)
filename = sprintf( 'Upper_Bound_weights_.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

%% Test chaning both upper & lower
rng(0);
pcer = 5e-2;
p = 1e3;
mu = -sort(abs(randn(p,1)),'ascend');
l = 1e3;
%Old:
% min_weights =((0:4:9)/10)';
% max_weights =1+((1:4:10)/10)';
min_weights = [0.25, 0.5];
max_weights =[1.5, 1.7];

w_m = zeros(length(min_weights)*length(max_weights),p);
for i=1:length(min_weights)
    for  j=1:length(max_weights)
        min_weight  = min_weights(i);
        max_weight  = max_weights(j);
        [w_m((i-1)*length(max_weights)+j,:)] = monotone_weights_sub(mu,pcer,min_weight,max_weight,l,1);
    end
end
rng(2);
a = {'-','--','-.',':'};
figure, hold on
h1 = plot(mu,w_m(1,:),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(mu,w_m(2,:),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
h3 = plot(mu,w_m(3,:),'linewidth',4,'color',rand(1,3));
set(h3,'LineStyle',a{3});
h4 = plot(mu,w_m(4,:),'linewidth',4,'color',rand(1,3));
set(h4,'LineStyle',a{4});
xlabel('Prior Mean');
ylabel('Weight');
set(gca,'fontsize',20)
xlim([min(mu),max(mu)]);
filename = sprintf( 'Upper_and_Lower_Bound_weights_.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf)
