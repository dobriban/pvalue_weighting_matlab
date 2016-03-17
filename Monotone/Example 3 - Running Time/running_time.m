%test monotone weights: Running time
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 3 - Running Time');
addpath '../../Code' '../../Code/Helper Code'

%% running time
K = 5;
num_exp = 50;
%% time Monotone weights
p_array = 10.^(2:K);
times = zeros(length(p_array),num_exp);
lb = zeros(length(p_array),num_exp);
rng(1);
for j=1:length(p_array)
    p = p_array(j);
    for i=1:num_exp
        fprintf('Experiment %d/%d for dimension %d.\n', i, num_exp,p)
        pcer = rand/10;
        min_weight = rand/10;
        %max_weight = 10+10*rand;
        mu = -abs(randn(p,1));
        tic
        [w_m,lb(j,i)] = monotone_weights_sub(mu,pcer,min_weight);
        times(j,i) = toc;
    end
end
%%
rng(2);
figure
%use delta method to scale standard errors to log scale
errorbar(log10(p_array),log10(mean(times')),2*sqrt(var(times'))./mean(times'),'linewidth',2,'color',rand(1,3))
xlabel('Log10(Dim)')
ylabel('Log10(Sec)')
xlim([min(log10(p_array)),max(log10(p_array))]);
set(gca,'fontsize',20)
filename = sprintf( 'Timing_small.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% check distibution of running times
hist(log(times(4,:)))
%% how many times break loop?
%9 out of 50
sum(lb')