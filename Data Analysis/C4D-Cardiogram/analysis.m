cd('C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\C4D-Cardiogram');
addpath '../../Code'
%%
prior = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011'; 
current = 'Cardiogram-CAD.Nat.Genet.2011';

%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis_flexible(prior,current,sensitivity_analysis);