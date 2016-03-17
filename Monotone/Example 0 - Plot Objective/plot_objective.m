%plot objective
cd('C:\Dropbox\Projects\Flexible p-value weighting\Monotone\Example 0 - Plot Objective');
addpath '../../Code' '../../Code/Helper Code'

%%
mu_array = [-1 -3];
q_array = [0.05,0.1];
a = {'-','--','-.',':'};

rng(2);
figure, hold on
for i=1:length(mu_array)
    for j=1:length(q_array)
        mu  = mu_array(i);
        q  = q_array(j);
        f   = @(w) normcdf(norminv(q.*w)-mu);
        [x,y]= fplot(f,[0,1/q]);
     h=plot(x,y,'linewidth',4,'color',rand(1,3));
     set(h,'LineStyle',a{(i-1)*2+j});
    end
end
xlabel('q');
ylabel('f');
legend('mu=-1 q=0.05', 'mu=-1 q=0.1', 'mu=-3 q=0.05', 'mu=-3 q=0.1', 'location','SouthEast');

set(gca,'fontsize',20)
filename = sprintf( 'Plot_Objective.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%saveTightFigure(gcf,'Plot_Objective.pdf')