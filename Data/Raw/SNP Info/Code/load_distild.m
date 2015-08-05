%Process the DistiLD data into MATLAB format
%keeps only the LD blocks
fileID = fopen('lds.tsv');
formatSpec = '%s %*s %*s*';

k = 0;
LD_chr = [];
LD_left_end = [];
LD_right_end = [];

%while k<10
while  ~feof(fileID)
  k = k+1;
  C = textscan(fileID,formatSpec,1,'Delimiter','\t');
  
  formatSpec1 = 'chr%d:%d-%d%*[^\n]';
  D = textscan(char(C{1}),formatSpec1);
  LD_chr = [LD_chr; D(1)];
  LD_left_end = [LD_left_end; D(2)];
  LD_right_end = [LD_right_end; D(3)];
end


fclose(fileID);

LD_chr = cell2mat(LD_chr);
LD_left_end = cell2mat(LD_left_end);
LD_right_end = cell2mat(LD_right_end);

%% delete trailing empty cells (possibly sex chromosomes)
L = length(LD_chr);
for i=1:L
  if isempty(LD_chr{i}) break
  end
end

LD_chr = cell2mat(LD_chr(1:i-1));
LD_left_end = cell2mat(LD_left_end(1:i-1));
LD_right_end = cell2mat(LD_right_end(1:i-1));

%%
save('lds.mat','LD_chr','LD_left_end','LD_right_end');

