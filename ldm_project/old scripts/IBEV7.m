%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
%  - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function out = IBEV7(data,error)

    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

% calculate the observed data
%     x = normrnd(data(:,1),error(:,1));
%     y = normrnd(data(:,2),error(:,2));
    x = data(:,1);
    y = data(:,2); 
    z = [x,y];
    
    [~, i1] = sort(z);
    [~, i2] = sort(i1);
    fobs = diff(i2,[],2)./n;
    
    
    % calculate expected data
    mu = nanmean(data,2);
    sigma = nansum(error,2)/2;
    g = normrnd(mu,sigma);
    h = normrnd(mu,sigma);
    
    j = [g,h];
    
    [~, i1] = sort(j);
    [~, i2] = sort(i1);
    fexp = diff(i2,[],2)./n;
    
    f = @(x) mean(abs(x - median(x)));
    out = f(fobs) - f(fexp);
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    