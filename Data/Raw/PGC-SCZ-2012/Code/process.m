%code to process data
%to be run in its own folder

%cd('C:\Dropbox\Weighted New\Bayes\Data Analysis\Data\Raw\PGC-SCZ-2012\Code\');
addpath '..'
fid = fopen('/../pgc.scz.full.2012-04.txt');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
%snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
formatSpec = 'rs%f %f %f %c %c %f %f %f %f %f %f';
num = size(Out,1)-1;
SNP = zeros(num,1);
chr = zeros(num,1);
bp = zeros(num,1);
A1 = char(zeros(num,1));
A2 = char(zeros(num,1));
or = zeros(num,1);
se = zeros(num,1);
P = zeros(num,1);
problem_snp = zeros(num,1);

tic
for i=2:num+1
    C = textscan(Out{i},formatSpec,'Delimiter','\t');
    if isempty(C{1})
        problem_snp(i-1) = 1;
    else
        SNP(i-1) = C{1};
        chr(i-1) = C{2};
        bp(i-1) = C{3};
        A1(i-1) = C{4};
        A2(i-1) = C{5};
        or(i-1) = C{6};
        se(i-1) = C{7};
        P(i-1) = C{8};
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
%%
N = 21856*ones(length(SNP),1);
%%
save('../../Processed/PGC-SCZ-2012.Nat.Genet.2012.mat','SNP','chr','bp','A1','A2','or','se','P','N');
