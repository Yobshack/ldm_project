%% Kyobi Skutt-Kakaria
% 01.20.2018
% Harvard University
% de Bivort Lab
% this function calculates the relative timing of turns compared to the light->dark transitions
% The function assigns every turn to a transition but doesnt sort based on whether a transition is
% light to dark or dark to light. the transition time vector can be used for that.


function [relTime, relTransType, ttOut] = calculateRelativeTime(tStamps,lightDat)

% generate transition time vector
tt = [0; cumsum(lightDat(1:(end-1),1))];

% subset by transitions that change light condition
tt0 = find([0; diff(lightDat(:,2))] ~= 0);
tl = lightDat(tt0,2);
tt = tt(tt0);
ttOut = lightDat(tt0,:);

% this cell array will contain vectors that will contain the relative time stamp for each turn, 
% negative means it comes before a transition and positive is after
relTime = zeros(length(tStamps),1);
% this cell will contain a vector that assigns the turn to either a light to dark transition or a
% dark to light transition: 1 = dark to light, 0 = light to dark
relTransType = relTime;

    for jj = 1:length(tStamps)
        tDiff = tStamps(jj) - tt;
        
        % time until next transition
        t1 = max(tDiff(tDiff < 0));
        % time since last transition
        t2 = min(tDiff(tDiff >= 0));
        if isempty(t2)
            relTime(jj) = tStamps(jj);
        elseif  (abs(t1) < t2)  
            relTime(jj) = t1;
            relTransType(jj) = tl(sum(tStamps(jj) >= tt)+1);
        else
            relTime(jj) = t2;
            relTransType(jj) = tl(sum(tStamps(jj) >= tt));
        end
        

    end
    
    relTransType = logical(relTransType);
    