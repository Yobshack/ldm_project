%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 01.30.2018 - Updated
% Harvard University
% de Bivort Lab
% this calculates all the behavioral metric for every fly and appends them to the cell array
% generated previously.

% 1.30.2018 added a filter to attempt to eliminate the transition periods

function [dataMat,errorMat] = generateBehavioralMetrics3(data,lightDat,times,bins,turnFilter)

% create choosable data parameters
tempTimes = times;
light = {'light','dark'};

% initialize arrays for store data and error
dataMat = nan(1,5,length(light),length(times),sum(bins>0));
errorMat = dataMat;

for kk = 1:sum(bins>0)
    for ii = 1:length(tempTimes)
        for jj = 1:length(light)
            idx1 = data{2} > tempTimes{ii}(1) & data{2} <= tempTimes{ii}(2);
            
            % normally i set these cumulative sums to [0; vect(1:(end-1))], should I have done
            % that here?
            accumDat = cumsum(lightDat(:,1));
            accT = [0; accumDat(1:(end-1))];
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
            
            if length(bins) > 2
                
                if jj == 1
                    intensity = nan(length(tStamps),1);
                    for gg = 1:length(tStamps)
                        intensity(gg) = lightDat(sum(tStamps(gg) > accT),3);
                    end
                    
                    
                    idx = intensity == bins(kk);
                    
                    turns = turns(idx);
                    walls = walls(idx);
                    tStamps = tStamps(idx);
                end
            end
            
                    this is important with fast transitions but less important for big bins
                    [tRel tRelType] = calculateRelativeTime(tStamps,lightDat);
            %         timeBin = sum(lightDat(accumDat>tempTimes{ii}(1) & accumDat<=tempTimes{ii}(2),1));
            
                    idx = (tRel < 0) & tRel > -600;%| (tRel > 2 & tRelType);
            
                    if range(data{2})/(tempTimes{ii}(2)-tempTimes{ii}(1)) < 0.5
                        continue
                    end
            
            idx = true(length(turns),1);
            
            if length(turns) > turnFilter
                
                nt = turnrate(turns(idx),tStamps(idx),lightDat);
                tb = turnbias(turns(idx));
                sw = switchiness(turns(idx),lightDat,tStamps(idx));
                wd = wallDist(walls(idx),turns(idx));
                cl = clumpiness(tStamps(idx),lightDat);
                
                dataMat(1,:,jj,ii,kk) = [nt(1),tb(1),sw(1),wd(1),cl(1)];
                errorMat(1,:,jj,ii,kk) = [nt(2),tb(2),sw(2),wd(2),cl(2)];
                
            end
            
        end
    end
    
end

