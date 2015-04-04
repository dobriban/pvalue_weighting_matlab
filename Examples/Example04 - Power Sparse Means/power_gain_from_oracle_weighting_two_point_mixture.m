%% power gain from p-value weighting
% distribution of means: mixture of point mass at 0 and another value
cd('C:\Dropbox\Weighted New\pvalue_weighting_matlab\Examples\Example04 - Power Sparse Means');
addpath '../../Code'
addpath '../../Code/Helper Code'
%%
q = 1e-3; %this leads to 4-fold improvements

%q = 5*1e-2;
power_unif = @(mu,pi_1) (1-pi_1)*q+ pi_1.*normcdf(norminv( q) - mu);


[X,Y] = meshgrid( -5:.1:0, 0:1e-3:0.4);
% 
% figure
% subplot(1,2,1);
 Z1 = power_unif(X,Y);
% contour(X,Y,Z1,'ShowText','on')
% title('power of uniform weighting');
% xlabel('size of large effect');
% ylabel('proportion of large effects');
% 
% subplot(1,2,2);
 Z2 = power_opt(X,Y,q);
% contour(X,Y,Z2,'ShowText','on')
% title('power of optimal weighting');
% xlabel('size of large effect');
% ylabel('proportion of large effects');


%%
figure
[C, H] = contour(X,Y,Z2./Z1,'ShowText','on');
set (H, 'LineWidth', 1);
%title('power of optimal / power of uniform');
xlabel('M');
ylabel('\pi_1');
set(gca,'LooseInset',get(gca,'TightInset'))
%%
saveTightFigure(gcf,'Power_Uniform_vs_Optimal_Weighting_two_Point_Mixture.pdf')
%% grayscale
colormap(gray)
saveTightFigure(gcf,'Power_Uniform_vs_Optimal_Weighting_two_Point_Mixture_grayscale.pdf')

