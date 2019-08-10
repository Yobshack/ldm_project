%% Kyobi Skutt-Kakaria
% 02.03.2018
% 02.04.2018 updated
% removed the attempt to correct for turns being across multiple bins as i realized that it doesnt
% help with this statistic. I can do it with switchiness though.
% Harvard University
% de Bivort Lab
% this function calculates the mutual information between current turn and next turn for an individual


function out = mutualinformation(turns)



% mutual information calculation, this seems right for now
% calculation of error
p1b = nan(100,1);
for ii = 1:100
    idx = randi(length(turns),length(turns),1);
    if ~all(turns(idx)) && any(turns(idx)) 
    temp = ami(turns(:,1),turns(:,2),0);
    p1b(ii) = temp(1);
    end
end


% mutual information

out = [nanmean(p1b) nanstd(p1b)];





