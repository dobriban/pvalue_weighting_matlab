
%% assess effect of changing pfer & sigma
M_alpha = 3;
step_size = 0.3;
sigma_array = (0.1:step_size:M_alpha)';
pfer_array =  (0.1:step_size:M_alpha)';

mu0 = sqrt(N_current./N_prior).*Z_prior;

num_sig_pfer = zeros(length(sigma_array), length(pfer_array));
which_sig_pfer = zeros(J,1); %how many times are hypotheses significant

for i=1:length(sigma_array)
    for j=1:length(pfer_array);
        sigma0 = sigma_array(i)*ones(J,1);
        pcer0 = pfer_array(j)/J;
        w_r0 = regularized_weights(mu0,sigma0,pcer0); %compute regularized weights
        P_wr0 = P_current./w_r0;
        [h0]=bonferroni(P_wr0,q,report);
        
        num_sig_pfer(i,j) = sum(h0);
        h0(isnan(h0))=0; %set the Nans to zero;
        which_sig_pfer = which_sig_pfer + h0;
    end
end

%% save a figure
figure
imagesc(sigma_array, pfer_array,num_sig_pfer)
colormap(flipud(pink));
colorbar
xlabel('sigma')
ylabel('pfer')
title('Sensitivity of number of rejections');
colorbar('location','eastoutside')
set(gca,'position',[0.1 0.1 0.7 0.8],'units','normalized')
filename = './Results/Robustness/Sensitivity_of_num_rej_pfer_sigma.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);

%%
colheadings = {'SNP','times sig.'};
rowheadings  = {};
rej_r0 = find(which_sig_pfer>=1);
snp_rej0 = snp_overlap(rej_r0);
num_rej0 = which_sig_pfer(rej_r0);

data = [snp_rej0 num_rej0];
wid = 16;
fms = {'.2d'};

colsep = ' & ';
rowending = ' \\';

filename = [analysis_ID ' - Sensitivity_of_of_num_rej_pfer_sigma_LaTeX'];
fileID = fopen(['./Results/Robustness/' filename '.txt'],'w');

displaytable(data,colheadings,wid,fms,rowheadings,fileID,colsep,rowending);

num_tests0 = length(pfer_array)*length(sigma_array);
fileIfD = fopen(['./Results/Robustness/' filename '.txt'],'at+');
fprintf(fileID,'\nResults of robustness analysis while changing alpha and PFER\n');
fprintf(fileID, 'Total Number of Attempts: %d\n',num_tests0);
fprintf(fileID,'The results displayed show the number of time each particular SNP was declared significant.\n');
fclose(fileID);

fprintf(['Saved Results to ./Results/Robustness/' filename '.txt\n']);