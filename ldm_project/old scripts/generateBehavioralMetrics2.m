%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 02.27.2018 - Updated
% Harvard University
% de Bivort Lab
% this calculates all the behavioral metric for every fly and appends them to the cell array
% generated previously.

% 1.30.2018 added a filter to attempt to eliminate the trans-transition time stamps for clumpiness
% and switchiness

function out = generateBehavioralMetrics2(data,lightDat)


% generates accumulated time throughout experiment
accumDat = cumsum(lightDat(:,1));
switch jj
    case 1
        idx2 = data{11};
        accumDat = accumDat(logical(lightDat(:,2)));
    case 2
        idx2 = ~data{11};
        accumDat = accumDat(logical(lightDat(:,2)));
    case 3
        idx2 = ones(length(idx1),1);
end
idx = idx1 & idx2;
turns = data{1}(idx);
walls = data{10}(idx);
tStamps = data{2}(idx);

% filter out turns that happen too closely together
diffs = diff(tStamps);
idx = [true; diffs<1];

turns = turns(~idx);
walls = walls(~idx);
tStamps = tStamps(~idx);

[tRel, tRelType] = calculateRelativeTime(tStamps,lightDat);
timeBin = sum(lightDat(accumDat>tempTimes{ii}(1) & accumDat<=tempTimes{ii}(2),1));

idx = (tRel > 5 & ~tRelType) | tRel < 0 | (tRel > 2 & tRelType);

if range(tStamps)/(tempTimes{ii}(2)-tempTimes{ii}(1)) < 0.5
    allDataCellArray{hh,jj,ii} = {[nan nan],[nan nan],[nan nan],[nan nan],[nan nan]};
    continue
end


if length(turns) > turnFilter
    
    nt = turnrate(turns(idx),timeBin,tStamps(idx),lightDat);
    tb = turnbias(turns(idx));
    sw = switchiness(turns(idx),lightDat,tStamps(idx));
    wd =  wallDist(walls(idx),turns(idx));
    cl = clumpiness(tStamps(idx),lightDat);
    out = {nt,tb,sw,wd,cl};
    
else
    out = {[nan nan],[nan nan],[nan nan],[nan nan],[nan nan]};
end


