%code to process SNP data
%to be run in its own folder

%cd('C:\Git\pvalue-weighting-gwas\Data\Raw\SNP Info\Code\');

%takes cca 3min on my desktop
addpath '..'
fid = fopen('/../snp141Common.txt');
tic
Out  = textscan(fid,'%s','delimiter',sprintf('\n'));
toc
fclose(fid);
Out = Out{1,1};
%%
%takes cca 20 min on my desktop
%need to generate 3 variables
%SNP chr pos

%format: 
%index chr# pos# pos#+1 rs# etc..
formatSpec = '%f chr%f %f  %f rs%f %s';
num = size(Out,1);
SNP = zeros(num,1);
chr= zeros(num,1);
pos= zeros(num,1);
problem_snp=zeros(num,1);

tic
for i=1:num
    C = textscan(Out{i},formatSpec,'Delimiter','\t');
    if isempty(C{1})||isempty(C{2})
        problem_snp(i) = 1;
    else
        SNP(i) = C{5};
        chr(i) = C{2};
        pos(i) = C{3};
    end
    K=1e4;
    if i/K==floor(i/K)
        str = sprintf('Processing SNP %d out of %d, time elapsed: %f sec.\n',i,num+1,toc);
        fprintf(str);
    end
end
toc

SNP = SNP(~problem_snp);
chr = chr(~problem_snp);
pos = pos(~problem_snp);
%%
save('../../Processed/snp141Common.mat','SNP','chr','pos');
