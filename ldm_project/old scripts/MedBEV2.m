%% Kyobi Skutt-Kakaria
% 03.03.2018 - Created
%  - Updated
% Harvard University
% de Bivort Lab
% this function is based off of IBEV3 but is used to calculate median effects instead of rank order

function out = MedBEV2(data,error)

    data = data(all(~isnan(data),2),:);
    error = error(all(~isnan(error),2),:);

    n = size(data,1);

% calculate the observed data
    x = normrnd(data(:,1),error(:,1));
    y = normrnd(data(:,2),error(:,2));
       
    z = [x,y];
%     
%     [~, i1] = sort(z);
%     [~, i2] = sort(i1);
%     fobs = diff(i2,[],2)./n;
%     
    % calculate median
%     fobs = diff(z,[],2);

     
    out = y; 
    
    %out2 = cat(3,fobs,fexp);
end
    
    
    