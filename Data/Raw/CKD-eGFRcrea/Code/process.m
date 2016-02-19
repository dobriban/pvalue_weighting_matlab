%code to process data
%to be run in its own folder

%cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Data\Raw\CKD-eGFRcrea\Code\');
addpath '..'
fid = fopen('/../CKDGen-eGFRcrea_meta_post.csv');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
%the format is
%'rsID,allele1,allele2,freqA1,direction,pval'
formatSpec = 'rs%f %c %c %s %c %f';
num = size(Out,1)-1;
SNP = zeros(num,1);
A1 = char(zeros(num,1));
A2 = char(zeros(num,1));
freqA1 = zeros(num,1);
direction =  char(zeros(num,1));
P = zeros(num,1);
problem_snp = zeros(num,1);

tic
for i=2:num+1
    C = textscan(Out{i},formatSpec,'Delimiter',',');
    if isempty(C{1})
        problem_snp(i-1) = 1;
    else
        SNP(i-1) = C{1};
        A1(i-1) = C{2};
        A2(i-1) = C{3};
        freqA1(i-1) = str2double(C{4});
        direction(i-1) = C{5};
        P(i-1) = C{6};
    end
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f.2 sec.\n',i,num+1,toc);
        fprintf(str);
    end
end
toc

SNP = SNP(~problem_snp);
A1 = A1(~problem_snp);
A2 = A2(~problem_snp);
freqA1 = freqA1(~problem_snp);
direction = direction(~problem_snp);
P = P(~problem_snp);

%%
N = 62237*ones(length(SNP),1);
%%
save('../../Processed/CKD-eGFRcrea.Kottgen.Nat.Genet.2010.mat','SNP','A1','A2','freqA1','direction','P','N');
