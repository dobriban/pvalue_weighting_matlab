%% Load Prior data
%loads overlapping SNPs with their chromosomal information
%as well as prior association information

%load snp info
load('../../Data/Processed/snp141Common.mat','SNP','chr','pos');

% clean data: find the overlapping SNPs; prior and current
[snp_1, ind1, ind2] = intersect(snp_overlap, SNP);
%subsample the other data files
chr = chr(ind2);
pos = pos(ind2);

num_snps_not_in_UCSC = length(snp_overlap) - length(snp_1);
% clear temporary
clear SNP

%End up with: snp_overlap, chr, pos, P_prior, N_prior
%% Cluster
clusters = cell(length(methods),1);
for i=1:length(methods)
    clusters{i} = distild_prune(snp_1,chr,pos,P_current(ind1),h(i,ind1));
end

%% Print to File
%print high level statistics
if ~(exist('./Results/','dir')==7)
mkdir('./Results/');
end

filename = [analysis_ID ' - clustered results'];
fileID = fopen(['./Results/' filename '.txt'],'w');

fprintf(fileID,['Results of P-value Weighting Analysis - Clustered Results\n\n']);
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
fprintf(fileID,['Number of Significant SNPs in Distinct LD Blocks: \n']);
for i=1:length(methods)
    str = sprintf('%s: %d\n',methods{i}, size(clusters{i},1));
    fprintf(fileID,str);
end
fprintf(fileID,'\n');

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
    for i=1:length(sigma_array)
        str = sprintf('P_thresh= %e\n',P_thresh(i));
        fprintf(fileID,str);
    end 
end
fprintf(fileID,'\n');

%print all SNPs with info
fprintf(fileID,['Lead Significant SNPs in each cluster: \n']);
for i=1:length(methods)
    str = sprintf('%s: %d SNPs\n',methods{i}, size(clusters{i},1));
    fprintf(fileID,str);
    
    %arrange SNPs 
    fprintf(fileID,'rsID\tchr\tposition\tcurrent_Pvalue\n');
    for j=1:size(clusters{i},1)
        fprintf(fileID,'rs%d\t%d\t%d\t%e\n',clusters{i}(j,1),clusters{i}(j,2),clusters{i}(j,3),clusters{i}(j,4));
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

filename = [analysis_ID ' - clustered results summary latex'];
fileID = fopen(['./Results/Latex/' filename '.txt'],'w');

fprintf(fileID,['Results of P-value Weighting Analysis -  clustered results \n\n']);
fprintf(fileID,['Prior Study: ' prior '\n']);
fprintf(fileID,['Current Study: ' current '\n']);
fprintf(fileID,'\n');


for i=1:length(methods)
    str = sprintf('%s\t',methods{i});
    fprintf(fileID,str);
end
fprintf(fileID,'\\\\ \n');

for i=1:length(methods)
    fprintf(fileID,'& \\multicolumn{1}{l|}{');
    fprintf(fileID,num2str(size(clusters{i},1)));
    fprintf(fileID,'}');
end
fprintf(fileID,'\\\\ \n');

fclose(fileID);
fprintf(['Saved Results to ./Results/' filename '.txt\n']);