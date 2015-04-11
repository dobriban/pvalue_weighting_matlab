%plot weights: Spjotvoll and regularized
cd('C:\Dropbox\Weighted New\pvalue_weighting_matlab\Examples\Example01 - Plot Weights/');
addpath '../../Code'
addpath '../../Code/Helper Code'
%% Plot Regularized weights
mini = -6;
maxi = 2;
pcer = 0.05;
min_sigma = 1e-1;
max_sigma = 4;
sigma_array = linspace(min_sigma, max_sigma,20);
epsi = 2*1e-1;
[X,Y] = meshgrid(mini:epsi:maxi, sigma_array);
Z = regularized_weights_unnorm(X,Y,pcer);
k = mean(mean(Z));
Z = 1/2*Z/k;

figure
colormap(brighten(jet,0.6))
surf(X,Y,Z)
xlim([mini,maxi])
ylim([min_sigma,max_sigma])
zlim([0,max(max(Z))])
xlabel('\eta');
ylabel('\sigma');
zlabel('w');
%%
saveTightFigure(gcf,'Regularized_Weights_Surface.pdf')
%% grayscale
colormap(brighten(gray,0.99))
%%
saveTightFigure(gcf,'Regularized_Weights_Surface_grayscale.pdf')

%% Contour

[X,Y] = meshgrid( mini:.01:maxi, 0.01:.01:4);
Z = regularized_weights_unnorm(X,Y,pcer);
k = mean(mean(Z));
Z = 1/2*Z/k;

figure
contour(X,Y,Z,'ShowText','on')
xlim([mini,maxi])
xlabel('\eta');
ylabel('\sigma');
colormap(brighten(jet,0.6))
%%
saveTightFigure(gcf,'Regularized_Weights_Contour.pdf')
%% grayscale
colormap(gray)
saveTightFigure(gcf,'Regularized_Weights_Contour_grayscale.pdf')

%% Not used: 
%% Plot Spjotvoll weights
c = 2;
pcer = 10/10^3;
w = @(mu) normcdf(c./mu+mu/2)/pcer;

figure
fplot(w,[-5,-0.01]);
title('Spjotvoll weights');
xlabel('mean');
ylabel('weight');

saveTightFigure(gcf,'Spjotvoll_Weights_Plot.pdf')

%% Plot Spjotvoll weights  - other modes
c_grid = 1;
mini  = -5;
maxi  = -0.01;
num = 1000;
grid = linspace(mini,maxi,num);
normw  = zeros(num,length(c_grid));

for i=1:length(c_grid)
    c = c_grid(i);
    pcer = 10/10^3;
    w = @(mu) normcdf(c./mu+mu/2)/pcer;
    
    w_grid = w(grid);
    normw(:,i) = w_grid/max(w_grid);
end

figure
hold on
plot(-grid,normw,'*');
%title('optimal weights');
xlabel('prior effect size');
ylabel('weight');
saveTightFigure(gcf,'Multi_Spjotvoll_Weights_Plot.pdf')
