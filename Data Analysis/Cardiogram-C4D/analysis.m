cd('C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\Cardiogram-C4D');
addpath '../../Code'
%%
prior = 'Cardiogram-CAD.Nat.Genet.2011';
current = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';

%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis_flexible(prior,current,sensitivity_analysis);