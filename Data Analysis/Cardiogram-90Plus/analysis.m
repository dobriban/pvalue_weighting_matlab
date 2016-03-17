cd('C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\Cardiogram-90Plus');
addpath '../../Code'
%%
prior = 'Cardiogram-CAD.Nat.Genet.2011';
current = '90Plus-Aging.Deelen.Hum Mol Gen.2014';

%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis_flexible(prior,current,sensitivity_analysis);