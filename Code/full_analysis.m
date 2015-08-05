function [h,clusters,snp_overlap,methods_used]= full_analysis(prior,current,sensitivity_analysis,methods,q,report,pcer,sigma_array,beta,P_thresh,...
    print_results,plot_results,analysis_ID, generate_tables)

addpath '../../Code/Weighting'
addpath '../../Code/External Helper Code'

[snp_overlap, P_prior, N_prior, Z_prior, P_current, N_current] = load_data(prior, current);
J = length(P_current);

if ~exist('methods','var')
    methods = {'Unweighted', 'Spjotvoll', 'Regularized', 'Exponential','Select top'};
end
methods_used = {};
if any(strcmp(methods,'Regularized'))
    if ~exist('sigma_array','var')
        sigma_array = [0.1 1 10];
    end
end
if any(strcmp(methods,'Exponential'))
    if ~exist('beta','var')
        beta = [1 2 4];
    end
end
if any(strcmp(methods,'Select top'))
    if ~exist('P_thresh','var')
        P_thresh = 10.^(-[2 4 6]);
    end
end
num_methods  = length(methods) + (length(sigma_array) - 1)  +...
    (length(P_thresh)-1);

%index of selected SNPs
h = zeros(num_methods,J);
current_method = 0;

if ~exist('q','var')
    q = 0.05;
end

if ~exist('report','var')
    report = 'yes';
end

%% (1) unweighted
if any(strcmp(methods,'Unweighted'))
    if strcmp(report,'yes')
        fprintf('Unweighted: ')
    end
    [h_u]=bonferroni(P_current,q,report);
    current_method = current_method + 1;
    h(current_method,:)=h_u;
    methods_used{1,length(methods_used)+1} =  'Unweighted';
end

%% (2) Spjotvoll weighting
if any(strcmp(methods,'Spjotvoll'))
    mu = sqrt(N_current./N_prior).*Z_prior;
    
    if ~exist('pcer','var')
        pcer = 1/J; %level used in weighting
    end
    
    w = spjotvoll_weights(mu,pcer);
    P_w = P_current./w;
    
    if strcmp(report,'yes')
        fprintf('Spjotvoll: ')
    end
    [h_s]=bonferroni(P_w,q,report);
    current_method = current_method + 1;
    h(current_method,:)=h_s;
    methods_used{1,length(methods_used)+1} =  'Spjotvoll';
end
%% (3) Regularized weighting
if any(strcmp(methods,'Regularized'))
    
    if ~exist('pcer','var')
        pcer = 1/J; %level used in weighting
    end
    sigma = sqrt(N_current./N_prior);
    w_r = zeros(J,length(sigma_array));
    for i=1:length(sigma_array)
        w_r(:,i) = regularized_weights(mu,sigma_array(i)*sigma,pcer);
        P_wr = P_current./w_r(:,i);
        if strcmp(report,'yes')
            str = sprintf('Regularized; sigma = %f: ',sigma_array(i));
            fprintf(str);
        end
        [h_r]=bonferroni(P_wr,q,report);
        current_method = current_method + 1;
        h(current_method,:)=h_r;
        methods_used{1,length(methods_used)+1} = str;
    end
end

%% (4) Exp weights
if any(strcmp(methods,'Exponential'))
    mu = sqrt(N_current./N_prior).*Z_prior;
    w_exp = zeros(J,length(beta));
    for i=1:length(beta)
        w_exp(:,i) = exp(beta(i)*abs(mu));
        w_exp(:,i) = J*w_exp(:,i)/sum(w_exp(:,i));
        P_wexp  = P_current./w_exp(:,i);
        if strcmp(report,'yes')
            str = sprintf('Exponential; beta = %f: ',beta(i));
            fprintf(str);
        end
        [h_wexp]=bonferroni(P_wexp,q,report);
        current_method = current_method + 1;
        h(current_method,:)=h_wexp;
        methods_used{1,length(methods_used)+1} = str;
    end
end

%% (5) Select according to p-value threshold, weight them uniformly
if any(strcmp(methods,'Select top'))
    for i=1:length(P_thresh)
        index = find(P_prior<P_thresh(i));
        if strcmp(report,'yes')
            str = sprintf('Select top; P_thresh = %e: ',P_thresh(i));
            fprintf(str);
        end
        [h_top]=bonferroni(P_current(index),q,report);
        current_method = current_method + 1;
        h(current_method,index)=h_top;
        methods_used{1,length(methods_used)+1} = str;
    end
end

%% update the methods
methods = methods_used;
%% post-analysis
clusters = [];
if ~exist('print_results','var')
    print_results = 1;
end
if ~exist('analysis_ID','var')
    analysis_ID = [prior ' - ' current];
end
if print_results==1
    run('post_analysis.m')
end

%% plotting
if ~exist('plot_results','var')
    plot_results = 0;
end
if plot_results==1
    run('post_analysis_plotting.m')
end

%% sensitivity analysis
if ~exist('sensitivity_analysis','var')
    sensitivity_analysis = 0;
end
if (sensitivity_analysis==1)
    run('robust.m')
end
