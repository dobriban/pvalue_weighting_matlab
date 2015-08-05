%% Sensitivity analysis for iGWAS

addpath '../../Code/Robustness Analysis'

%% plot: scatter z-scores; QQ-Plot or bivariate normality
run('scatter_z_scores.m');

%% test Spjotvoll weights computed by using a two-parameter family of 
%%piecewise linear shrinkers on the means
run('sensitivity_shrinkers.m');

%% assess effect of changing pfer & sigma on regularized weights
run('sensitivity_alpha_pfer.m');

%% Assess sensitivity of results to parameter choices
%vary prior variance when computing regularizedc weights
run('sensitivity_prior_variance.m');

%% assess effect of correlations on modelling
%assume the prior and current means are correlated, and vary correlation
%parameter
run('Correlated_Effects_Weights.m');
