%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
% 03.03.2018 - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function out = IBEV10(data,error)


    %f = @(x) x./nanmean(x);
    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

% calculate the observed data
    %x = normrnd(data,error);
    x = data;
    z = [x(:,1)-x(:,2), x(:,3)-x(:,4)];
    
%     idx = abs(z(:,1)) > nanmean(error(:,1));
%     z(~idx,1) = nan;
%     idx = abs(z(:,2)) > nanmean(error(:,2));
%     z(~idx,2) = nan;
    % calculate expected data
    mu = [nanmean(data(:,1:2),2),nanmean(data(:,3:4),2)];
    sigma = [nanmean(error(:,1:2),2),nanmean(error(:,3:4),2)];
    g = normrnd(mu,sigma);
    h = normrnd(mu,sigma);
    
    j = g-mu;
    
    %f = @(x) nanmean(abs(x-median(x)));
    
    f = @(x) mad(x);
    temp = f(z) - f(j);
    
    out = (temp(2)/(temp(1)+0.01))*temp(2);
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    