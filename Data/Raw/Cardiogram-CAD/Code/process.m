%code to process data
%to be run in its own folder

%cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Data\Raw\Cardiogram-CAD\Code\');
addpath '..'
fid = fopen('/../CARDIoGRAM_GWAS_RESULTS.TXT');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
formatSpec = 'rs%f chr%f:%f %c %c %f %f %f %f %f %f %f %s';
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
        N(i-1) = C{11} + C{12};
        log_odds(i-1) = C{9};
        log_odds_se(i-1)=C{10};
        P(i-1) = C{7};
    end
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f.2 sec.\n',i,num+1,toc);
        fprintf(str);
    end
end
toc


%%
save('../../Processed/Cardiogram-CAD.Nat.Genet.2011.mat','SNP','chr','pos','A1','A2','freqA1','log_odds','log_odds_se','P','N');
