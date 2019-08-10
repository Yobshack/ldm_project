%% Kyobi Skutt-Kakaria
% 02.10.2018 - Created
%  - Updated
% Harvard University
% de Bivort Lab
% this is the function fed into a bootstrapping function to calculate the change in behavioral
% metrics it is a monte carlo metric that calculates the average difference between rank position of
% the same individual across conditions

function out = IBEV8(data,error)

    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

% calculate the observed data
    x = normrnd(data,error);
       
    z1 = x(:,1:2);
    z2 = x(:,3:4);
    
    [~, i1] = sort(z1);
    [~, i2] = sort(i1);
    
    [~, i3] = sort(z2);
    [~, i4] = sort(i3);
    
    
    fobs = [diff(i2,[],2),diff(i4,[],2)]./n;
    
    
    % calculate expected data
    x = normrnd(data,error);
    y = normrnd(data,error);
    
    z1 = [x(:,1), y(:,1)];
    z2 = [x(:,3), y(:,3)];
    
    [~, i1] = sort(z1);
    [~, i2] = sort(i1);
    
    [~, i3] = sort(z2);
    [~, i4] = sort(i3);

    fexp = [diff(i2,[],2),diff(i4,[],2)]./n;
    
    f = @(x) mean(abs(x - median(x)));
    %f = @(x) mad(x);
    temp = f(fobs) - f(fexp);
    out = [temp(2)/temp(1) temp];
    %out2 = cat(3,fobs,fexp);
end
    
    
    