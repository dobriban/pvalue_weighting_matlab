function [h]=bonferroni(pvals,q,report)

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        %error('Some p-values are greater than 1.');
    end
end

if nargin<2,
    q=.05;
end

if nargin<3,
    report='no';
end

s=length(pvals);
thresh = q/s;
h = (pvals <= thresh);


if strcmpi(report,'yes'),
    n_sig=sum(pvals<=h);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a family-wise error rate of %f.\n',s,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a  family-wise error rate  of %f.\n',s,n_sig,q);
    end
end




