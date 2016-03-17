%this file outputs the results of the iGWAS
%compare the outputs of all the methods:
%Note: if you add new methods, need to edit only the section where the
%parameters are printed

%load snp info
load('../../Data/Processed/snp141Common.mat','SNP','chr','pos');
[~, ind1, ~] = intersect(snp_overlap, SNP);
num_snps_not_in_UCSC = length(snp_overlap) - length(ind1);


%print high level statistics
if ~(exist('./Results/','dir')==7)
    mkdir('./Results/');
end

filename = [analysis_ID ' - full results'];
fileID = fopen(['./Results/' filename '.txt'],'w');

fprintf(fileID,['Results of P-value Weighting Analysis\n\n']);
fprintf(fileID,['Prior Study: ' prior '\n']);
fprintf(fileID,['Current Study: ' current '\n']);
fprintf(fileID,'\n');

str = sprintf('Number of overlapping SNPs: %d\n',J);
fprintf(fileID,str);
str = sprintf('Number of overlapping SNPs not found in UCSC database: %d\n',num_snps_not_in_UCSC);
fprintf(fileID,str);
fprintf(fileID,'\n');

fprintf(fileID,['Weighted Bonferroni Methods Tested: \n']);
for i=1:length(methods)
    fprintf(fileID,[methods{i} '; ']);
end
fprintf(fileID,'\n');
str = sprintf('Family-wise Error Rate controlled at: %e\n',q);
fprintf(fileID,str);
fprintf(fileID,'\n');
fprintf(fileID,['Number of Significant SNPs: \n']);
for i=1:length(methods)
    str = sprintf('%s: %d\n',methods{i}, sum(h(i,:)));
    fprintf(fileID,str);
end
fprintf(fileID,'\n');

%parameter settings
fprintf(fileID,'Parameter settings:\n');
if exist('q','var')
    str = sprintf('q = %e\n',q);
    fprintf(fileID,str);
end
if exist('pcer','var')
    str = sprintf('pcer = %e\n',pcer);
    fprintf(fileID,str);
end
if exist('sigma_array','var')
    for i=1:length(sigma_array)
        str = sprintf('sigma = %e\n',sigma_array(i));
        fprintf(fileID,str);
    end
end


if exist('beta','var')
    str = sprintf('beta= %e\n',beta);
    fprintf(fileID,str);
end
if exist('P_thresh','var')
    for i=1:length(P_thresh)
        str = sprintf('P_thresh= %e\n',P_thresh(i));
        fprintf(fileID,str);
    end
end
if exist('lower_bounds','var')
    for i=1:length(lower_bounds)
        str = sprintf('lower_bounds= %e\n',lower_bounds(i));
        fprintf(fileID,str);
    end
end
if exist('upper_bounds','var')
    for i=1:length(upper_bounds)
        str = sprintf('upper_bounds= %e\n',upper_bounds(i));
        fprintf(fileID,str);
    end
end
fprintf(fileID,'\n');

%print all SNPs with info
fprintf(fileID,['Significant SNPs: \n']);
for i=1:length(methods)
    str = sprintf('%s: %d SNPs\n',methods{i}, sum(h(i,:)));
    fprintf(fileID,str);
    
    %arrange SNPs
    [snp_r,chr_r,pos_r] = snp_arrange(h(i,:),snp_overlap,SNP,chr,pos);
    snps_discovered = snp_overlap(h(i,:)==1);
    
    fprintf(fileID,'rsID\tchr\tposition\tcurrent_Pvalue\n');
    for j=1:length(snp_r)
        ind = find(snp_overlap ==snp_r(j));
        pval_r  = P_current(ind);
        fprintf(fileID,'rs%d\t%d\t%d\t%e\n',snp_r(j),chr_r(j),pos_r(j),pval_r);
    end
    fprintf(fileID,'\n');
end
fclose(fileID);
fprintf(['Saved Results to ./Results/' filename '.txt\n']);

%% Latex
%print high level statistics
if ~(exist('./Results/Latex','dir')==7)
    mkdir('./Results/Latex');
end

filename = [analysis_ID ' - full results summary latex'];
fileID = fopen(['./Results/Latex/' filename '.txt'],'w');

fprintf(fileID,['Results of P-value Weighting Analysis\n\n']);
fprintf(fileID,['Prior Study: ' prior '\n']);
fprintf(fileID,['Current Study: ' current '\n']);
fprintf(fileID,'\n');


for i=1:length(methods)
    str = sprintf('%s\t',methods{i});
    fprintf(fileID,str);
end
fprintf(fileID,'\\\\ \n');

%& \multicolumn{1}{l|}{number}
for i=1:length(methods)
    fprintf(fileID,'& \\multicolumn{1}{l|}{');
    fprintf(fileID,num2str(sum(h(i,:))));
    fprintf(fileID,'}');
end
fprintf(fileID,'\\\\ \n');



fclose(fileID);
fprintf(['Saved Results to ./Results/' filename '.txt\n']);
