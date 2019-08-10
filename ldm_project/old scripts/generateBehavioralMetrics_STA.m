%% Kyobi Skutt-Kakaria
% 01.31.2018 - Created
% 01.30.2018 - Updated
% Harvard University
% de Bivort Lab
% this calculates all the behavioral metric for every fly and generates transition triggered averages

function out = generateBehavioralMetrics(allDataCell,lightDat)

% create choosable data parameters
temps = {'low','high'};
tempTimes = {[0 80]*60,[120 200]*60};
light = {'light','dark'};
turnFilter = 30;

% Generate cell array with behavioral metric for each individual to be indexed later
allDataCellArray = cell(1,length(temps),length(light));
for hh = 1:size(allDataCell,1)
    data = allDataCell(hh,:);
    
    
    hh
    for ii = 1:length(temps)
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
            
            [tRel tRelType] = calculateRelativeTime(tStamps,lightDat);
            timeBin = sum(lightDat(accumDat>tempTimes{ii}(1) & accumDat<=tempTimes{ii}(2),1));
            
            idx = tRel > 10 | tRel < -5;
            
            turns = turns;
            walls = walls;
            tStamps = tStamps;
            timeVect = [-20:5:20 30];
            
                
            for kk = 1:(length(timeVect)-1)
                   
                    idx = tRel >= timeVect(kk) & tRel <= timeVect(kk+1);
                    
                    
                if sum(idx) > turnFilter
                    allDataCellArray{hh,jj,ii,kk} = {numturns(turns(idx)),turnbias(turns(idx)),...
                        switchiness(turns(idx)), wallDist(walls(idx),turns(idx)),...
                        clumpiness(tStamps(idx),lightDat)};
                
                else
     
                    allDataCellArray{hh,jj,ii,kk} = {nan,[nan nan],[nan nan],[nan nan],[nan nan]};
                end
            end
            
        end
    end
    
end

out = allDataCellArray;
