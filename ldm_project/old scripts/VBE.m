%% Kyobi Skutt-Kakaria
% 02.01.2018 - Created
% 02.02.2018 - Updated
% Harvard University
% de Bivort Lab

function out = VBE(data,error)

% calculate the observed data
    x = normrnd(data,error);
    % calculate expected data
    mu = nanmean(data);
    sigma = nansum(error,2)/2;
    g = normrnd(mu,sigma);
    
    out = mad(x,2) - mad(g,2);
end
    
    
    