%% Kyobi Skutt-Kakaria
% 01.26.2018 - Created
% 01.30.2018 - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between two individuals
% normalized by 1 + the average distance

function out = IBE(data,error)

% calculate the observed data
    x = normrnd(data(:,1),error(:,1));
    y = normrnd(data(:,2),error(:,2));
    
    f = @(x) (x-nanmean(x))./nanstd(x);
    x = f(x);
    y = f(y);
    
    fobs = abs(x - y);
    gobs = nanmean(x - y);
    aobs = 1;%0.5*nanmean(x + y);
    hobs = nanmean(fobs./(aobs+gobs));
    
    % calculate expected data
    mu = nanmean(data,2);
    sigma = sum(error)/2;
    g = normrnd(mu,sigma);
    h = normrnd(mu,sigma);
    g = f(g);
    h = f(h);
    
    fexp = abs(g - h);
    gexp = nanmean(g - h);
    aexp = 1;%0.5*nanmean(x + y);
    hexp = nanmean(fexp./(aexp+gexp));
    
    out = hobs - hexp;
end
    
    
    