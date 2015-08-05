% display jointly the snps selected in varying shrinkers and alpha/pfer

colheadings = {'SNP','shrinkers', 'PFER/sigma'};
rowheadings  = {};
rej_any = find(max((which_sig>=1),(which_sig_pfer>=1)));
snp_rej = snp_overlap(rej_any);
num_rej_shrink = which_sig(rej_any);
num_rej_pfer = which_sig_pfer(rej_any);

l = length(lower_knot_array);
num_tests_shrink = l*(l+1)/2;
num_tests_pfer = length(pfer_array)*length(sigma_array);

data = [snp_rej num_rej_shrink/num_tests_shrink num_rej_pfer/num_tests_pfer];
wid = 16;
fms = {'.2d','.2f','.2f'};

colsep = ' & ';
rowending = ' \\';

fileID = 1;
displaytable(data,colheadings,wid,fms,rowheadings,fileID,colsep,rowending);
