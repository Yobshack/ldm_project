%% Kyobi Skutt-Kakaria
% 01.20.2018
% Harvard University
% de Bivort Lab
% this function calculates the switchiness score for an individual
% CONSIDER CHANGING THIS TO MUTUAL INFORMATION, WHAT WOULD BE DIFFERENT? %
% Right now, this function calculates the probability of seeing a right turn following a right turn
% divided by the probability of right turns over all, a score of 1 means that you observe a right
% following a right with the same probability of seeing a right overall and argues against chaining
% of turns.

function out = switchiness(turns,lightDat,tStamps)

turnfilter = 40;

% generate transition time vector
tt = [0; cumsum(lightDat(1:(end-1),1))];

% subset by transitions that change light condition
tt0 = find([0; diff(lightDat(:,2))] ~= 0);
tt = [0; tt(tt0)];

% return turns that happen on either side of a time transition
ti = zeros(length(tStamps)-1,2);
for ii = 1:(length(tStamps)-1)
    ti(ii,:) = [sum(tStamps(ii) > tt),sum(tStamps(ii+1) > tt)];
end   


% number of right turns following right turns
rrturns = turns(1:(end-1)) + turns(2:end)==1;

% remove turns found across transitions;
rrturns(ti(:,1) ~= ti(:,2)) = [];
turns([true;ti(:,1) ~= ti(:,2)]) = [];

if size(rrturns,1) > turnfilter

    % make sums of r turn vector
    for ii = 1:100
        idx = randi(length(rrturns),length(rrturns),1);
        n1b(ii) = sum(rrturns(idx));
        n2b(ii) = nanmean(turns(idx));
        if n1b(ii) == 0
            n1b(ii) = nan;
            n2b(ii) = nan;
        end
    end

probR = n2b;
    
% normalization factor, need to figure out why this is the way it is
p2 = 2.*probR.*(1-probR).*length(turns);

% stat % inverting this toggles between switchiness and streakiness
% stat = n1b./p2;
stat = p2./n1b;

% switchiness score = ratio of right turns followed by right turns divided by overall R bias

out = [mean(stat) std(stat)];
else
    out = [nan nan];
end





