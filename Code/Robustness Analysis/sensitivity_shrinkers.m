% test a two-parameter family of piecewise linear shrinkers

%% obtain the rejected hypotheses for each setting
M = 10;
step_size = 1;
lower_knot_array = (0:step_size:M)';
upper_knot_array =  (0:step_size:M)';

% J = length(Z_prior);
% pcer = 1/J; %level used in weighting
% q = 0.05; report = 'yes';

num_sig = nan(length(lower_knot_array), length(upper_knot_array));
which_sig = zeros(J,1); %how many times are hypotheses significant
eps = 1e-3;

for i=1:length(lower_knot_array)
    for j=i:length(upper_knot_array);
        a = lower_knot_array(i); b = upper_knot_array(j);
        mu0 = piecewise_linear_shrinker(-Z_prior,a,b); %need to supply negative means
        mu0 = -sqrt(N_current./N_prior).*mu0;
        [w_r0,~,err] = spjotvoll_weights(mu0,pcer);
        if err==1
            break
        end
        P_wr0 = P_current./w_r0;
        [h0]=bonferroni(P_wr0,q,report);
        num_sig(i,j) = sum(h0);
        h0(isnan(h0))=0; %set the Nans to zero;
        which_sig = which_sig + h0;
    end
    if err==1
        break
    end
end

%% plot number of rejections

figure
imagesc(lower_knot_array, upper_knot_array,num_sig)
colormap(flipud(pink));
%colorbar
xlabel('Upper Knot')
ylabel('Lower Knot')
%title('Sensitivity of number of rejections');
%colorbar('location','eastoutside')
%set(gca,'position',[0.1 0.1 0.8 0.9],'units','normalized')
filename = './Results/Robustness/Sensitivity_of_num_rej_shrinkers.pdf';
saveTightFigure(gcf,filename);
fprintf(['Saved Results to ' filename '\n']);
%% write to file how many times each is rejected
%% in Latex format
if ~exist('generate_tables','var')
    generate_tables = 1;
end
if generate_tables == 1
    colheadings = {'SNP','times sig.'};
    rowheadings  = {};
    rej_r0 = find(which_sig>=1);
    snp_rej0 = snp_overlap(rej_r0);
    num_rej0 = which_sig(rej_r0);
    
    data = [snp_rej0 num_rej0];
    wid = 16;
    fms = {'.2d'};
    
    colsep = ' & ';
    rowending = ' \\';
    
    filename = ['Sensitivity_of_num_rej_shrinkers_LaTeX'];
    fileID = fopen(['./Results/Robustness/' filename '.txt'],'w');
    
    displaytable(data,colheadings,wid,fms,rowheadings,fileID,colsep,rowending);
    
    l = length(lower_knot_array);
    num_tests = l*(l+1)/2;
    
    fileIfD = fopen(['./Results/Robustness/' filename '.txt'],'at+');
    fprintf(fileID,'\nResults of robustness analysis while using piecewise linear shrinkers\n');
    fprintf(fileID, 'Total Number of Parameter settings tested: %d\n',num_tests);
    fprintf(fileID,'The results displayed show the number of time each particular SNP was declared significant.\n');
    fclose(fileID);
    fprintf(['Saved Results to ./Results/Robustness/' filename '.txt\n']);
end

