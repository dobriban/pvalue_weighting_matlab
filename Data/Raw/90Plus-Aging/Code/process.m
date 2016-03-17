%code to process data
%to be run in its own folder
%% load file
cd('C:\Dropbox\Projects\Flexible p-value weighting\Data\Raw\90Plus-Aging\Code\');
addpath '..'
fid = fopen('/../ddu139supp_data_file.txt');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
%% create some basic vars
Out = Out{1,1};
num = size(Out,1)-1;
K = 1e5;
%%
%SNP	Chr_hg18build36	Position_hg18build36	Non_Effect_Allele	Effect_Allele	Direction_of_effect_85	Pvalue_85	Direction_of_effect_90	Pvalue_90
%rs10	7	92221824	A	C	-	.742382675247881	+	.195299931370681
formatSpec = 'rs%f %f %f %c %c %c %f %c %f';
SNP = zeros(num,1);
chr =  zeros(num,1);
pos =  zeros(num,1);
A1 = char(zeros(num,1));
A2 = char(zeros(num,1));
dir_85 = char(zeros(num,1));
P_85 =  zeros(num,1);
dir = char(zeros(num,1));
P =  zeros(num,1);
problem_snp = zeros(num,1);

tic
for i=2:num+1
    C = textscan(Out{i},formatSpec,'Delimiter','\t');
    if isempty(C{1})||isempty(C{9})
        problem_snp(i-1) = 1;
    else
        SNP(i-1) = C{1};
        chr(i-1) =  C{2};
        pos(i-1) =  C{3};
        A1(i-1) = C{4};
        A2(i-1) = C{5};
        dir_85(i-1) = C{6};
        P_85(i-1) =  C{7};
        dir(i-1) = C{8};
        P(i-1) =  C{9};
    end
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f.2 sec.\n',i,num+1,toc);
        fprintf(str);
    end
end
toc

%% clear problem SNPs
sum(problem_snp) %10394

SNP = SNP(~problem_snp);
chr =  chr(~problem_snp);
pos = pos(~problem_snp);
A1 = A1(~problem_snp);
A2 = A2(~problem_snp);
dir_85 = dir_85(~problem_snp);
P_85 =  P_85(~problem_snp);
dir = dir(~problem_snp);
P = P(~problem_snp);

%% sample size from the paper
N =  (5406+16121)*ones(length(SNP),1);
%%
tic
save('../../../Processed/90Plus-Aging.Deelen.Hum Mol Gen.2014.mat','SNP','chr','pos','A1','A2','dir_85','P_85','dir','P','N');
toc