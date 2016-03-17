cd('C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\Lipids-SCZ');
addpath '../../Code'
%%
prior = 'Lipids-TG-ONE.Teslovich.Nature.2010';
current = 'PGC-SCZ-2012.Nat.Genet.2012';

%% 
sensitivity_analysis = 0;
[h,clusters,snp,methods] = full_analysis_flexible(prior,current,sensitivity_analysis);