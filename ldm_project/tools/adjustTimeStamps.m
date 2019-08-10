function out = adjustTimeStamps(flyTracks)

for kk = 1:size(flyTracks.centroid,3)
   % This finds a distance from the center of the maze
   
    x = flyTracks.centroid(:,:,kk);
    x = x - min(x);
   
    % this rotates all the mazes between 65 and 120 so that the analysis works on those
    % mazes as well
    if kk > 64
        x(:,2) = max(x(:,2)) - x(:,2);
        x(:,1) = max(x(:,1)) - x(:,1);
    end
    
    % return maximum yValue
    maxes = max(x);
    xmax = maxes(1);
    ymax = maxes(2);
    
    % This finds a distance from the center of the maze
    distFromCent = sqrt(sum((x - [xmax*0.5 ymax*0.667]).^2,2));
    meanShape = mean(x(distFromCent < 10,:));
    distFromCent = sqrt(sum((x - meanShape).^2,2));
   
    % this calculates a change in distance from the center at every frame
    deltaDist = [0; diff(distFromCent)];

    fTTurns = flyTracks.rightTurns(:,kk);
    allTInd = find(~isnan(fTTurns));
    
    newInd = nan(length(allTInd),1);
    for jj = 1:length(allTInd)
    
    backThresh = 15;
    distBack = 0;
    count = 1;
    while distBack < backThresh
        if count == (allTInd(jj)-1)
            distBack = 15;
        else
        distBack = distBack + abs(deltaDist(allTInd(jj)-count));
        count = count+1;
        end
    end
    

    [~,b] = min(distFromCent((allTInd(jj)-count):allTInd(jj)));
    
    newInd(jj) = allTInd(jj) - (count - b + 1);
    end
       
    flyTracks.rightTurns(newInd,kk) = flyTracks.rightTurns(allTInd,kk);
    flyTracks.rightTurns(allTInd,kk) = nan;
        

end


out = flyTracks;