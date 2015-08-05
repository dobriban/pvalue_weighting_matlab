
%% Assess sensitivity of results to parameter choices
%vary prior variance when computing regularizedc weights
M = 10;
step_size = 1;
sigma2_array = (0.1:step_size:M)'; %array of variances

num_sig = zeros(length(sigma2_array),1);
    
for i=1:length(sigma2_array)
    s2 =  sigma2_array(i);
    
    mu0 = sqrt(N_prior).*s2./(1+N_prior.*s2).*Z_prior;
    mu0 = sqrt(N_current).*mu0;
    
    sigma0 = s2/(1+N_prior*s2)*ones(J,1);
    sigma0 = sqrt(N_current.*sigma0);

    w_r0 = regularized_weights(mu0,sigma0,pcer); %compute regularized weights
    P_wr0 = P_current./w_r0;
    [h0]=bonferroni(P_wr0,q,report);
    num_sig(i) = sum(h0);
end
%%
%  plot(mu0, w_r0, '*')
%%
figure
plot(sigma2_array,num_sig);
xlabel('Prior Variance')
ylabel('Number of Rejections')
filename = './Results/Robustness/Sensitivity_of_num_rej_prior_variance.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);