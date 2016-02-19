%code to process data
%to be run in its own folder

%cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Data\Raw\Lipids-TG-ONE\Code\');
addpath '..'
fid = fopen('/../TG_ONE_Eur.tbl.sorted.txt');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
formatSpec = 'rs%f %c %c %f %f %f %c %s';
num = size(Out,1)-1;
SNP = zeros(num,1);
A1 = char(zeros(num,1));
A2 = char(zeros(num,1));
Weight = zeros(num,1);
GC_Zscore = zeros(num,1);
P = zeros(num,1);
Overall = char(zeros(num,1));
problem_snp = zeros(num,1);

tic
for i=2:num+1
    C = textscan(Out{i},formatSpec,'Delimiter','\t');
    if isempty(C{1})
        problem_snp(i-1) = 1;
    else
        SNP(i-1) = C{1};
        A1(i-1) = C{2};
        A2(i-1) = C{3};
        Weight(i-1) = C{4};
        GC_Zscore(i-1) = C{5};
        P(i-1) = C{6};
        Overall(i-1) = C{7};
    end
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f.2 sec.\n',i,num+1,toc);
        fprintf(str);
    end
end
toc

% SNP = SNP(~problem_snp);
% A1 = A1(~problem_snp);
% A2 = A2(~problem_snp);
% freqA1 = freqA1(~problem_snp);
% direction = direction(~problem_snp);
% P = P(~problem_snp);
sum(problem_snp)
%% from the paper
N =  96598*ones(length(SNP),1);
%%
save('../../Processed/Lipids-TG-ONE.Teslovich.Nature.2010.mat','SNP','A1','A2','Weight','GC_Zscore','P','Overall','N');
