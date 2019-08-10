%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 01.30.2018 - Updated
% Harvard University
% de Bivort Lab
% this calculates all the behavioral metric for every fly and appends them to the cell array
% generated previously.

% 1.30.2018 added a filter to attempt to eliminate the transition periods

function out = generateBehavioralMetrics(allDataCell,lightDat,times)

% create choosable data parameters
tempTimes = times;
% tempTimes = {[0 240*60] [0 240*60]};
light = {'light','dark'};
turnFilter = 100;


% Generate cell array with behavioral metric for each individual to be indexed later
allDataCellArray = cell(1,length(tempTimes),length(light));
for hh = 1:size(allDataCell,1)
    data = allDataCell(hh,:);


hh
for ii = 1:length(tempTimes)
    uInt = unique(lightDat(:,3));
    for jj = 1:length(light)
        idx1 = data{2} > tempTimes{ii}(1) & data{2} <= tempTimes{ii}(2);
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
        idx = [true; diffs<2];
        
        turns = turns(~idx);
        walls = walls(~idx);
        tStamps = tStamps(~idx);
        
        [tRel tRelType] = calculateRelativeTime(tStamps,lightDat);
        timeBin = sum(lightDat(accumDat>tempTimes{ii}(1) & accumDat<=tempTimes{ii}(2),1));
        
        idx = (tRel > 5 & ~tRelType) | tRel < 0 | (tRel > 2 & tRelType);
        
        if range(tStamps)/(tempTimes{ii}(2)-tempTimes{ii}(1)) < 0.5
            allDataCellArray{hh,jj,ii} = {[nan nan],[nan nan],[nan nan],[nan nan],[nan nan]};
            continue
        end
       
        turns = turns;
        walls = walls;
        tStamps = tStamps;
        
        if length(turns) > turnFilter

%             nt = turnrate(turns(idx),timeBin);
%             tb = turnbias(turns(idx));
%             sw = switchiness(turns(idx));
%             wd = wallDist(walls(idx),turns(idx));
%             cl = clumpiness(tStamps(idx),lightDat);
        allDataCellArray{hh,jj,ii} = {turnrate(turns(idx),timeBin,tStamps(idx),lightDat),...
            turnbias(turns(idx)),...
            switchiness(turns(idx),lightDat,tStamps(idx)),...
            wallDist(walls(idx),turns(idx)),clumpiness(tStamps(idx),lightDat)};
        
        else
            allDataCellArray{hh,jj,ii} = {[nan nan],[nan nan],[nan nan],[nan nan],[nan nan]};
        end
        
    end
end

end

out = allDataCellArray;
