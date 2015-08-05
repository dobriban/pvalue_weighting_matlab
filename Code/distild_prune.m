function clusters = distild_prune(snp,chr,pos,P,h)
%Cluster SNPs based on Distild data
%Inputs
%snp - list of snps
%chr, pos - genomic information
%P - p-calues, used to pick best snps
%h - indicatior set of which snps to work with (say L total)

%Outputs
%clusters - a 4xL array, containing the genomic information for the
%selected snps

% %Testing code
% snp = snp_overlap;
% P = P_current;
% h = h_u;

%start by applying subsetting:
snp = snp(h==1);
P = P(h==1);
pos = pos(h==1);
chr = chr(h==1);


%% Load distiLD LD data
load('../../Data/Processed/DistiLD.lds.mat','LD_chr','LD_left_end','LD_right_end');

%%prune
%% sort by chromosome
[chr, ind] = sort(chr);
P = P(ind);
pos = pos(ind);
snp = snp(ind);
%% walk through chromosomes
C = max(chr);

pruned_snp = [];
pruned_P = [];
pruned_chr = [];
pruned_bp = [];

for i=1:C
  %get the current data on chromosome i
  ind = find(chr==i);
  current_snps = snp(ind);
  current_pos = pos(ind);
  current_p = P(ind);
  
  %sort based on bp
  [current_pos, ind1] = sort(current_pos);
  current_snps = current_snps(ind1);
  current_p = current_p(ind1);
  
  %find LD blocks on current chromosome
  ind_LD = find(LD_chr==i);
  %start the left index
  snp_index_l = 1;
  
  %step through the LD blocks
  for j=1:length(ind_LD)
    %find the basebair endpoints of the LD block
    left = LD_left_end(ind_LD(j));
    right = LD_right_end(ind_LD(j));
    
    %find the first snp inside the LD block (to pass left endpoint)
    while(snp_index_l<length(current_snps)) && (current_pos(snp_index_l)<left)
      snp_index_l = snp_index_l+1;
    end
    
    %find the last snp inside LD block
    snp_index_r = snp_index_l;
    while (snp_index_r<=length(current_snps)) && (current_pos(snp_index_r)<=right)
      snp_index_r = snp_index_r+1;
    end
    
    %if there is at least one GWAS snp in the LD block
    %find the best one
    if (snp_index_l<snp_index_r)
      snp_index_r = snp_index_r-1;
      if (current_pos(snp_index_l)>=left) && (current_pos(snp_index_r)<=right)
       current_LD_p = current_p(snp_index_l:snp_index_r);
        current_LD_snp = current_snps(snp_index_l:snp_index_r);
        current_LD_pos = current_pos(snp_index_l:snp_index_r);
        
        best_snp = find(current_LD_p == min(current_LD_p));
        ind = i*ones(length(best_snp),1);
        
        pruned_snp = [pruned_snp; current_LD_snp(best_snp)];
        pruned_P = [pruned_P; current_LD_p(best_snp)];
        pruned_chr = [pruned_chr; ind];
        pruned_bp = [pruned_bp; current_LD_pos(best_snp)];
      end
    end
    %update left pointer
    snp_index_l = snp_index_r;
    
  end
end 

clusters = [cast(pruned_snp,'double'), cast(pruned_chr,'double'),cast(pruned_bp,'double'), pruned_P];
