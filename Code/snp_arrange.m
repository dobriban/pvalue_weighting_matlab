function [snp_r,chr_r,pos_r] = snp_arrange(h_u,snp_overlap,SNP,chr,pos)
%arrange SNPs in order of chromosomes and location in genome;
%this script is used in the pipeline of comparing the results of methods

%Inputs
%h_u - index of significant SNPs in snp_overlap
%snp_current,chr,pos - annotation information

%Outputs
%snp_r,chr_r,pos_r - sorted results

rej = find(h_u==1);
%store the (SNP,Chromosome,Position) combo
snp_r = zeros(length(rej),1);
chr_r = zeros(length(rej),1);
pos_r = zeros(length(rej),1);
for i=1:length(rej)
    snp_r(i) = snp_overlap(rej(i));
    ind = find(SNP==snp_r(i));
    if ~isempty(ind)
        chr_r(i) = chr(ind);
        pos_r(i) = pos(ind);
    else
        chr_r(i) = NaN;
        pos_r(i) = NaN;
    end
end

%sort by chromosome
[chr_r, ind] = sort(chr_r);
snp_r = snp_r(ind);
pos_r = pos_r(ind);

%for each chromosome, sort snps by genome location
unique_chr = unique(chr_r);
for i=1:length(unique_chr);
    this_chr_ind = find(chr_r==unique_chr(i));
    [pos_r_new, ind] = sort(pos_r(this_chr_ind));
    pos_r(this_chr_ind) = pos_r_new;
    snp_sub = snp_r(this_chr_ind);
    snp_r(this_chr_ind) = snp_sub(ind);
    chr_sub = chr_r(this_chr_ind);
    chr_r(this_chr_ind) = chr_sub(ind);
end