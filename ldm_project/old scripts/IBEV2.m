%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
%  - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function out = IBEV2(data,error)

% calculate the observed data
    x = normrnd(data(:,1),error(:,1));
    y = normrnd(data(:,2),error(:,2));
    
    f = @(x) (x-nanmean(x))./nanstd(x);
    x = f(x);
    y = f(y);
    
    fobs = x - y;
    
    % calculate expected data
    mu = nanmean(data,2);
    sigma = nansum(error,2)/2;
    g = normrnd(mu,sigma);
    h = normrnd(mu,sigma);
    g = f(g);
    h = f(h);
    
    fexp = g - h;
    
    out = mad(fobs,2) - mad(fexp,2);
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    