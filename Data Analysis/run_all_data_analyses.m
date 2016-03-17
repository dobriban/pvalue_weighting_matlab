%Run all data analyses
dir = 'C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\';

analyses = {
'C4D-90Plus',
'C4D-Cardiogram',
'C4D-eGFRcrea',
'Cardiogram-90Plus',
'Cardiogram-C4D',
'eGFRcrea-C4D',
'Lipids-SCZ',
'SCZ-90Plus',
'SCZ-Lipids'
'Lipids-90Plus',
'eGFRcrea-90Plus'
};
%%
t = tic;
for i=1:length(analyses)
    a = analyses(i); a = a{:};
    cd(['C:\Dropbox\Projects\Flexible p-value weighting\Data Analysis\' a]);
    run('analysis.m');
end
toc(t);
