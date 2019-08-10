%% Kyobi Skutt-Kakaria
% 02.01.2018 - Created
% 02.02.2018 - Updated
% Harvard University
% de Bivort Lab

function [out,out2] = VBE2(data,error)

% calculate the observed data
    x = normrnd(data,error);
    % calculate expected data
    mu = repmat(nanmean(data),size(x,1),1);
    sigma = error;
    g = normrnd(mu,sigma);
    
    temp = mad(x,1) - mad(g,1);
    out = [temp(3:4)./temp(1:2) temp];
end
    
    
    