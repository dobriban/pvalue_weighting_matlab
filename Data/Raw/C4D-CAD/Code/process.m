%code to process data
%to be run in its own folder

%cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Data\Raw\C4D-CAD\Code\');

addpath '..'
fid = fopen('/../C4D_CAD_DISCOVERY_METAANALYSIS_UPDATE.TXT');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
%the format is
%SNP	CHR_POSB36	REFERENCE_ALLELE	OTHER_ALLELE	REFERENCE_ALLELE_FREQ	N_CASE	N_CONTROL	LOG_ODDS	LOG_ODDS_SE	PVALUE	HET_PVALUE	MODEL
%tab-separated
formatSpec = 'rs%f chr%f_%f %c %c %f %f %f %f %f %f %f %s';
num = size(Out,1)-1;
SNP = zeros(num,1);
A1 = char(zeros(num,1));
A2 = char(zeros(num,1));
freqA1 = zeros(num,1);
P = zeros(num,1);
problem_snp = zeros(num,1);
chr= zeros(num,1);
pos= zeros(num,1);
N= zeros(num,1);
log_odds= zeros(num,1);
log_odds_se= zeros(num,1);

tic
for i=2:num+1
    C = textscan(Out{i},formatSpec,'Delimiter','\t');
    if isempty(C{1})
        problem_snp(i-1) = 1;
    else
        SNP(i-1) = C{1};
        chr(i-1) = C{2};
        pos(i-1) = C{3};
        A1(i-1) = C{4};
        A2(i-1) = C{5};
        freqA1(i-1) = C{6};
        N(i-1) = C{7} + C{8};
        log_odds(i-1) = C{9};
        log_odds_se(i-1)=C{10};
        P(i-1) = C{11};
    end
    K=1e4;
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f.2 sec.\n',i,num+1,toc);
        fprintf(str);
    end    
end
toc


%%
save('../../Processed/C4D-CAD.C4D.Consortium.Nat.Genet.2011.mat','SNP','chr','pos','A1','A2','freqA1','log_odds','log_odds_se','P','N');
