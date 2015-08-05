cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Analysis\Cardiogram-C4D');
addpath '../../Code'
%%
prior = 'Lipids-TG-ONE.Teslovich.Nature.2010';
current = 'PGC-SCZ-2012.Nat.Genet.2012';
%
fig = stratified_qq_plot(prior,current); 
%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis(prior,current,sensitivity_analysis);