function [snp_overlap, P_prior, N_prior, Z_prior, P_current, N_current]= load_data(prior, current)
% loads uniformly processed data for p-value weighting analysis
% both prior and current data sets must be stored in a specific folder
%
% and they must have P,SNP,N variables

load(['../../Data/Processed/' prior '.mat'],'P','SNP','N'); 
P_prior = P;
snp_prior = SNP;
N_prior = N;

load(['../../Data/Processed/' current '.mat'],'P','SNP','N'); 
P_current = P;
snp_current = SNP;
N_current = N;

% clean data: find the overlapping SNPs; prior and current
[snp_overlap, ind1, ind2] = intersect(snp_prior, snp_current);
%subsample the other data files
P_prior = P_prior(ind1);
N_prior = N_prior(ind1);
P_current = P_current(ind2);
N_current = N_current(ind2);

Z_prior = norminv(P_prior/2);
%get rid of zero Z-s
nonzero_z = Z_prior<0;
P_prior = P_prior(nonzero_z);
Z_prior = Z_prior(nonzero_z);
N_prior = N_prior(nonzero_z);
P_current = P_current(nonzero_z);
N_current = N_current(nonzero_z);
snp_overlap = snp_overlap(nonzero_z);
