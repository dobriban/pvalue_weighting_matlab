cd('C:\Git\pvalue_weighting_matlab\Examples\Example09 - Optimization Plots');
addpath '../../Code'
addpath '../../Code/Helper Code'
%% plot summand in objective function
eta_a = [-0.1 -1];
sigma_a = [1 2];
q = 1/20;

a = {'-','--',':','-.'};
figure, hold on

col = colormap(gray(10)); %old version: col = colormap(gray(5));
%col = colormap(jet);

for i=1:length(eta_a)
  for j=1:length(sigma_a)
    eta = eta_a(i);
    sigma  = sigma_a(j);
    gamma = sqrt(1+sigma^2);   
    generic_term = @(w) normcdf( (norminv(q*w)-eta)/gamma);
    h = ezplot(generic_term,[0,1/q]);
    ylim([0,1]);
    set(h,'LineWidth',1);  %# Sets the line width to 1
    %set(h,'LineWidth',2); for talk
    c = a(2*(i-1)+j);
    set(h,'LineStyle',c{1}); 
    set(h,'color',col(2*(i-1)+j,:));
  end
end
%h_legend = legend('\eta=-0.1,\sigma=1', '\eta=-0.1,\sigma=2', '\eta=-1,\sigma=1', '\eta=-1,\sigma=2', ...
  %'location', 'southeast');
%set(h_legend,'FontSize',7);
title([]);
set(gca,'LooseInset',get(gca,'TightInset'))

%% plot
saveTightFigure(gcf,'objective_function.pdf')

%% plot
saveTightFigure(gcf,'objective_function_talk.pdf')


%% grayscale
saveTightFigure(gcf,'objective_function_grayscale.pdf')


%% Not used: 
%%
addpath '../External Helper Code'
%plot generic term - as a function of c
eta = -1;
sigma = 1;
gamma = sqrt(1+sigma^2);
lambda = 1.2;
generic_term = @(c) normcdf( (c-eta)/gamma) - lambda*normcdf(c);

h = ezplot(generic_term,[-10,10]);
set(h,'LineWidth',2);  %# Sets the line width to 2

%%
filename = sprintf( './generic_term.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);