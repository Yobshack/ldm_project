%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
% 03.03.2018 - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function [out,out2] = IBEV13(data,error)


    %f = @(x) x./nanmean(x);
    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

    f = @(x) (data-nanmean(data,1))./nanstd(data,[],1);
% %     
    error = error./nanstd(data,[],1);
    data = f(data);

    
% calculate the observed data
    %x = normrnd(data,error);
    x = data;
    z = [corr(x(:,1),x(:,2),'type','Pearson'), corr(x(:,3),x(:,4))];
    
%     idx = abs(z(:,1)) > nanmean(error(:,1));
%     z(~idx,1) = nan;
%     idx = abs(z(:,2)) > nanmean(error(:,2));
%     z(~idx,2) = nan;
%     % calculate expected data
%     mu = [nanmean(data(:,1:2),2),nanmean(data(:,3:4),2)];
%     sigma = [nanmean(error(:,1:2),2),nanmean(error(:,3:4),2)];
%     g = normrnd(mu,sigma);
%     h = normrnd(0,sigma);
%     
%     x = abs(z);
%     j = [(2*mean(error(:,1:2),2))./sqrt(pi),(2*mean(error(:,3:4),2))./sqrt(pi)];
%     
%     %f = @(x) nanmean(abs(x-median(x)));
% 
%     f = @(x) var(x);
%     temp = mad(z) ;%- nanmean(j);%-nanmedian(j);
    
%     out = (temp(2).^2/(temp(1)+0.01));
    out = sqrt(z(2).^2/z(1).^2);
    out2 = n;
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    