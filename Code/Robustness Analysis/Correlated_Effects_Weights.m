%% assess effect of correlations on modelling
%assume the prior and current means are correlated, and vary correlation
%parameter
%%
L = 4;
rho_array = linspace(0.9,1,L);

num_rejected = zeros(L,1);

figure


for i=1:L
  rho = rho_array(i);
  mu0 = rho*sqrt(N_current./N_prior).*Z_prior;
  sigma0 = sqrt(N_current./N_prior + N_current*(1-rho^2));
  w_r0 = regularized_weights(mu0,sigma0,pcer); %compute regularized weights
  P_wr0 = P_current./w_r0;
  
  [h_r0]=bonferroni(P_wr0,q,report);
  num_rejected(i)= sum(h_r0);
  subplot(2,2,i)
  plot(mu0, w_r0, '*')  ; hold on;
  xlabel('mean'); ylabel('weight');
    str = sprintf('rho = %f', rho);
  title(str);
end

%%
% num_selected = sum(w_r0>1);
% weight_of_selected = sum(w_r0(w_r0>1))/J;
% fprintf('Number Selected: %d. Weight of Selected: %f. \n',num_selected, weight_of_selected);

%%
filename = './Results/Robustness/Correlated_Effects_Weights.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);