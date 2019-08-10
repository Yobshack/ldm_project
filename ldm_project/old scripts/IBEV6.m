%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
% 03.13.2018 - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function out = IBEV6(data,error)

    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

% calculate the observed data
    x = normrnd(data(:,1),error(:,1));
    y = normrnd(data(:,2),error(:,2));
       
    z = x-y;
    
    % calculate expected data
    mu = nanmean(data,2);
    sigma = nanmean(error,2);
    g = normrnd(mu,sigma);
    h = normrnd(mu,sigma);
    
    j = g-h;
    
    f = @(x) x;
    out = f(z) - f(j);
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    