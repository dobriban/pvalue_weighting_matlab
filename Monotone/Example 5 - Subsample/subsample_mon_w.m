% compare monotone weights with the subsampling method
%vary the size of the problem (number of means) and the size of the
%subsample; record results
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 5 - Subsample');
addpath '../../Code' '../../Code/Helper Code'

%% Note: version 2 is used in the paper
%% Version 1:  keep a fixed number of means
savefigs=1; a = {'-','--','-.',':'};
pcer = 5e-3;
J_arr = [1e2, 1e3, 1e4, 2*1e4];%10.^(3:1:7);
l_arr = J_arr;
err = zeros(length(J_arr),length(l_arr));
min_weight = 0;
max_weight = Inf;
version  = 1;

for i=1:length(J_arr)
    J = J_arr(i);
    for j=1:i-1
        rng(0);  mu = -abs(randn(J,1));
        l = l_arr(j);
        
        tic
        w_m = monotone_weights(mu,pcer,min_weight,max_weight);
        toc
        tic
        w_m_s = monotone_weights_sub(mu,pcer,min_weight,max_weight,version,l);
        toc
        err(i,j) = mean(abs(w_m-w_m_s));
        
        if savefigs==1
            rng(2);
            log_err = log10(abs(w_m-w_m_s));
            [mu1,ind] = sort(mu);
            figure, hold on
            h1 = plot(mu1, log_err(ind),'linewidth',4,'color',rand(1,3));
            set(h1,'LineStyle',a{1});
            legend([h1],{'Err'},'location','Best')
            xlabel('mu')
            ylabel('log_{10} Err');
            set(gca,'fontsize',20)
            filename = sprintf( './Img/diff_mon_subsampling_J_=_%d_l_=_%d.png',J,l);
            saveas(gcf, filename,'png');
            fprintf(['Saved Results to ' filename '\n']);
            close(gcf)
        end
    end
end

%% Version 2: keep all means above a fixed gap
savefigs=1; a = {'-','--','-.',':'};
pcer = 5e-3;
J_arr = [1e2, 1e3, 1e4, 2*1e4];%10.^(3:1:7);
%J_arr = 1e2;
err = zeros(length(J_arr),1);
min_weight = 0;
max_weight = Inf;
version  = 2;
for i=1:length(J_arr)
    J = J_arr(i);
    rng(0);  mu = -abs(randn(J,1));
    tic
    w_m = monotone_weights(mu,pcer,min_weight,max_weight);
    toc
    tic
    w_m_s = monotone_weights_sub(mu,pcer,min_weight,max_weight,version);
    toc
    err(i) = mean(abs(w_m-w_m_s));
    
    if savefigs==1
        %relative error
        rng(2);
        log_err = log10(max(abs(w_m-w_m_s)./abs(w_m), 10^(-300)));
        [mu1,ind] = sort(mu);
        figure, hold on
        h1 = plot(mu1, log_err(ind),'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        legend([h1],{'Err'},'location','Best')
        xlabel('mu')
        xlim([min(mu1), max(mu1)]);
        ylabel('log_{10} Err');
        set(gca,'fontsize',20)
        filename = sprintf( './Img/diff_mon_subsampling_J_=_%d_version_2_rel.png',J);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
        
        %absolute error
        rng(2);
        log_err = log10(max(abs(w_m-w_m_s), 10^(-300)));
        [mu1,ind] = sort(mu);
        figure, hold on
        h1 = plot(mu1, log_err(ind),'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        legend([h1],{'Err'},'location','Best')
        xlabel('mu')
        xlim([min(mu1), max(mu1)]);
        ylabel('log_{10} Err');
        set(gca,'fontsize',20)
        filename = sprintf( './Img/diff_mon_subsampling_J_=_%d_version_2_abs.png',J);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
end
