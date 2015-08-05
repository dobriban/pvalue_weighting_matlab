if ~(exist('./Results/Robustness/','dir')==7)
mkdir('./Results/Robustness/');
end


%% create scatter plot of effect sizes 
Z_current = norminv(P_current/2);
x = [Z_current,Z_prior];

%% density estimate
[bandwidth,density,X,Y]=kde2d(x); 
% plot the data and the density estimate 
figure
contour3(X,Y,density,50), hold on 
plot(x(:,1),x(:,2),'r.','MarkerSize',5) 
title('Scatter and Density of Z-scores');
xlabel('Current Z-scores'); ylabel('Prior Z-scores');
filename = './Results/Robustness/Scatter_Density_z_Scores.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);

%% assess bivariate normality
r =Z_prior.^2+Z_current.^2;
theta = atan(Z_current./Z_prior);
J = length(Z_prior);
unif = linspace(0,1,J);
expo = - 2 * log ( 1 - unif);
%%
figure
subplot(1,2,1)
qqplot(expo,r)
xlabel('Chi Squared'); ylabel('Z-score ratio')
ylim([0,max(r)]);

subplot(1,2,2)
qqplot(linspace(0,pi/2,J), theta)
xlabel('Uniform'); ylabel('Angle')
ylim([0,max(theta)]);
figtitle('QQplots of Z-score ratio and Angle');

filename = './Results/Robustness/QQPlot_r_theta.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);
