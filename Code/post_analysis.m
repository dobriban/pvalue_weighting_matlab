%% Post-analysis for iGWAS


%% write to file the significant snps
% including genomic information
run('generate_output_with_locus_info.m');

%% DistiLD prune & write to file the significant snps
run('distild_pruning_results.m');

