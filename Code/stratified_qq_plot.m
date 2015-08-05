function fig = stratified_qq_plot(prior,current,func,funcname, lower_thresh, upper_thresh)
%%
addpath '../../Code/Helper Code'

[~, P_prior,~, ~, P_current, ~] = load_data(prior, current);

if ~exist('func','var')
    func  = @(x) - log10(x);
    funcname = '-log10';
end
if ~exist('lower_thresh','var')
    lower_thresh = [0 1 2 3];
end
if ~exist('upper_thresh','var')
    upper_thresh = -log10(5*10^(-8));
end

f_prior = func(P_prior);
f_current = func(P_current);

%%
rng(0);
m = Inf;
M = -Inf;
fig = figure; hold on
f=@(x) x;

h =zeros(3,length(lower_thresh));
for i=1:length(lower_thresh)
    ind = (f_prior>lower_thresh(i));
    theoretical_quantiles = (1:length(ind))/(length(ind)+1);
    theoretical_quantiles = func(theoretical_quantiles);
    m = min(min(theoretical_quantiles),m); M = max(max(theoretical_quantiles),M);
    h(:,i) = qqplot_modified(theoretical_quantiles,f_current(ind));
    col = rand(1,3);
    set(h(1,i),'markeredgecolor',col,'marker','.','markersize',12);
end
x = ezplot(f,[m,M]);
set(x, 'Color', 'k');
ylabel([funcname '(P) ' prior]);
xlabel('Expected quantiles');
title([ current ' | ' prior ]);
ylim([min(f_current),upper_thresh]);
xlim([m,M]);

legendCell = {};
for i=1:length(lower_thresh)
    legendCell{1,length(legendCell)+1} = num2str(lower_thresh(i), [funcname '(P_{prior})>=%-d']);
end
legend(h(1,:),legendCell,'location','best');


%%
if ~(exist('./Results/','dir')==7)
    mkdir('./Results/');
end

filename = sprintf( './Results/Stratified_QQ_Plot %s.png', funcname);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);