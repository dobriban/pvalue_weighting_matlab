cd('C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\eGFRcrea-C4D');
addpath '../../Code'
%%
prior = 'CKD-eGFRcrea.Kottgen.Nat.Genet.2010';
current = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';
%%
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis_flexible(prior,current,sensitivity_analysis);