cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Analysis\Cardiogram-C4D');
addpath '../../Code'
%%
prior = 'Cardiogram-CAD.Nat.Genet.2011';
current = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';
%
fig = stratified_qq_plot(prior,current); 
%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis(prior,current,sensitivity_analysis);