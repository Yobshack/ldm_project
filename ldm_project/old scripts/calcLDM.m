%% Kyobi Skutt-Kakaria
% 01.26.2018 - Created
% 01.26.2018 - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between two individuals
% normalized by 1 + the average distance

function out = calcLDM(data,error,ii)

    x = normrnd(data,error);

    fobs = abs(x(:,1) - x(:,2));
    gobs = nanmean(x(:,1) - x(:,2));
    aobs = (nanmean(nanmean(x,2),1));
    hobs = nansum(fobs./(aobs+gobs))*(1/size(x,1));
    
    x = [normrnd(nanmean([data(:,1), data(:,2)],2), sqrt(error(:,1).^2 + error(:,2).^2)/2),...
        normrnd(nanmean([data(:,1), data(:,2)],2), sqrt(error(:,1).^2 + error(:,2).^2)/2)];
    fexp = abs(x(:,1) - x(:,2));
    gexp = nanmean(x(:,1) - x(:,2));
    aexp = (nanmean(nanmean(x,2),1));
    hexp = nansum(fexp)./(aexp+gexp)*(1/size(x,1));
    
    out = hobs - hexp;
end
    
    
    