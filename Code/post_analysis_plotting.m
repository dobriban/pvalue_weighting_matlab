if ~(exist('./Results/','dir')==7)
mkdir('./Results/');
end

%% p-values
n = floor(sqrt(length(unique(P_prior))));
figure,
hist(P_prior,n)
str = ['prior p-values: ' prior];
title(str)
filename = sprintf( './Results/Pvalues_prior.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf);

%%
n = floor(sqrt(length(unique(P_current))));
figure,
hist(P_current,n)
str = ['current p-values: ' current];
title(str)
filename = sprintf( './Results/Current_Pvalues.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf);

%%  scatter z-scores
Z_current = norminv(P_current/2);
figure,
scatter(Z_prior,Z_current,'.')
title('Scatter of Z-scores');
xlabel(['Prior Z: ' prior]); ylabel(['Current Z: ' current]);
filename = sprintf( './Results/Scatter_zScores.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
close(gcf);

%% plot the weights
if any(strcmp(methods,'Regularized'))
    i = min(2,size(w_r,2));
    figure, plot(mu,w_r(:,i),'*')
    xlabel('mean'); ylabel('weight');
    str = fprintf('Regularized; sigma = %f', sigma_array(i));
    filename = './Results/Regularized_Weights_on_Data.pdf';
    saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
end
close(gcf);

%% min weight & actual weights
%only do this if have spjotvoll, regularized & exp
if any(strcmp(methods,'Spjotvoll'))&&any(strcmp(methods,'Exponential'))&&any(strcmp(methods,'Regularized'))
        i = min(2,size(w_r,2));
    acceptable_fwer = 0.05;
    min_p = acceptable_fwer/J;
    min_weight = P_current/min_p;
    feasible = (min_weight<1e4);
    figure, hold on
    scatter(Z_prior(feasible),log10(min_weight(feasible)))
    str  = sprintf('Minimum weight required for rejection, FWER=%f',acceptable_fwer);
    title(str);
    xlabel('Prior Z-scores'); ylabel('log(weight)');
    
    % plot the weights
    
    non_negligible = (w>0.1);
    plot(Z_prior(non_negligible),log10(w(non_negligible)),'.r')
    
    non_negligible = (w_r(:,i)>0.1);
    plot(Z_prior(non_negligible),log10(w_r(non_negligible,i)),'.g')
    
    non_negligible = (w_exp>0.1);
    plot(Z_prior(non_negligible),log10(w_exp(non_negligible)),'.k')
    
    legend('MinWeight', 'Spjotvoll','Regularized','Exponential','Location','SouthWest');
    filename = sprintf( './Results/Min_weight_required_and_actual_weights.png');
    saveas(gcf, filename,'png');
    
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf);
end
