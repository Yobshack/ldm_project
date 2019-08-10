function [out] = calculateWallDistV2(turnDirs,flyTracks)


for ii = 1:size(flyTracks.centroid,3)
    x = flyTracks.centroid(:,:,ii);
    x = x - min(x);
    

    
    
    % this rotates all the mazes between 65 and 120 so that the analysis works on those
    % mazes as well
    if ii > 64
        x(:,2) = max(x(:,2)) - x(:,2);
        x(:,1) = max(x(:,1)) - x(:,1);
    end
    
    % return maximum yValue
    maxes = max(x);
    xmax = maxes(1);
    ymax = maxes(2);
    
    % return shape previously calculated
    xshape =  x(~isnan(x(:,1)),:);
    shape = xshape(boundary(xshape,1),:);
    deviations = std(x);
    
    if isempty(shape) || any(deviations < 8.5) ||  any(deviations > 10)
        continue
    end
    shape = [smooth(interp(shape(:,1),5),50) smooth(interp(shape(:,2),5),50)];
    
    % This finds a distance from the center of the maze
    distFromCent = sqrt(sum((x - [xmax*0.5 ymax*0.667]).^2,2));
    meanShape = mean(x(distFromCent < 10,:));
    distFromCent = sqrt(sum((x - meanShape).^2,2));
    
    % this calculates a change in distance from the center at every frame
    deltaDist = [0; diff(distFromCent)];
    
    % this sets the "approach" thresholds
    farThresh = max(distFromCent)*0.667;
    closeThresh = max(distFromCent)*0.15;
    
    % assigns frames as "approach" frames if they are within a set distance of the choice point and
    % the fly is moving towards the choice point
    approaches = (distFromCent < farThresh & distFromCent > closeThresh) & (sign(deltaDist) == -1);
    
    % the meat of this script the idea is to calculate which wall the animal is closer to in every
    % approach frame
    dir = nan(length(turnDirs.tSequence{ii}),1);
    
    for jj = 1:length(turnDirs.tSequence{ii})
        % extracts time stamp for each turn calculated from calculateTurnDirsV2.m
        tTim = turnDirs.tTime{ii}(jj);
        
        % finds frames that are approach frames and also fall within a certain amount of time before
        % the fly passes through the choice point
        [~,b1] = min(abs(flyTracks.tStamps - tTim));
        distBack = zeros(1,1);
        count = 1;
        while distBack < farThresh
            distBack = distBack + abs(deltaDist(b1-count));
            count = count+1;
        end
        b2 = b1 - count;
        cent = x(b2:b1,:);
        distCent = distFromCent(b2:b1);
        distCentApp = distCent(any(diff(cent) > 0.1,2));
        % this adds a movement threshold that says that the animal needs to have moved 0.1 pixel for
        % that frame to count, this is just to avoid pauses getting overweighted in the position
        % distribution
        centApp = cent(any(diff(cent) > 0.1,2),:);
        % finds the average change in both x and y coordinate, this will be used to determine
        % which arm the animal is in
        if size(centApp,1) > 2
            deltaMean = mean(diff(centApp));
        elseif size(centApp,1) <= 1
            continue
        else
            deltaMean = diff(centApp);
        end
        
        avgMean = mean(centApp);

        
        % for each centroid, this finds a point on the boundary outline that is closest to the
        % centroid at that frame
        
        [a,b] = min(sqrt((centApp(:,1)-shape(:,1)').^2 + (centApp(:,2)-shape(:,2)').^2));
        distances = [centApp(:,1) - shape(b,1), centApp(:,2) - shape(b,2)];
        s = sign(centApp - meanShape);
        armMat(:,1) = s(:,1) == -1 & s(:,2) == 1;
        armMat(:,2) = s(:,1) == 1 & s(:,2) == 1;
        armMat(:,3) = s(:,2) == -1;
        
        wallDist = nan(size(centApp,1),2);
        for ll = 1:size(centApp,1)
            dist = centApp(ll,:) - shape;
            if arm == 1
                leftSideDist = dist(dist(:,2) < 0,:);
                rightSideDist = dist(dist(:,2) >= 0,:);
            elseif arm == 2
                leftSideDist = dist(dist(:,1) < 0,:);
                rightSideDist = dist(dist(:,1) >= 0,:);
            else
                leftSideDist = dist(dist(:,2) < 0,:);
                rightSideDist = dist(dist(:,2) >= 0,:);
            end
            wallDist(ll,:) = [min(sqrt(sum(leftSideDist.^2,2))) min(sqrt(sum(rightSideDist.^2,2)))];
        end
        
        centDistances{jj} = distCentApp;
        wallDistances{jj} = wallDist;
        turnDirection{jj} = repmat(turnDirs.tSequence{ii}(jj),size(centApp,1),1);
        
        % this figure is to show that the algorithm is working as intended and marking approaches as
        % either on the left or right side
        if dir(jj) == 1
            c = 'blue';
        elseif dir(jj) == 0
            c = 'red';
        elseif isnan(dir(jj))
            c = 'black';
        end
        %         scatter(centApp(:,1),centApp(:,2),'MarkerFaceColor',c,'Marker','.','MarkerEdgeColor',c,...
        %             'MarkerEdgeAlpha',0.5)
        %         hold on
    end

    
    allCentDistances = cat(1,centDistances{:});
    % This finds a distance from the center of the maze
    allWallDistances = cat(1,wallDistances{:});
    allTurnDirection = cat(1,turnDirection{:});

    
% assigns to the turnDirs struct
turnDirs.tWallDist{ii} = dir;

idx = allCentDistances > 3 & allTurnDirection & allCentDistances < 10;

[N,~,~,binX,binY] = histcounts2(allCentDistances(idx), allWallDistances(idx,2));
tic
histMat = nan(size(N))
for mm = 1:size(N,1)
    toc
    for nn = 1:size(N,2)
        idx = binX == mm & binY == nn;
        histMat(mm,nn) = nansum(allTurnDirection(idx))/length(allTurnDirection(idx));
    end
end

imagesc(histMat)
shg


scatter(allCentDistances(idx),allWallDistances(idx,2))


end

out = turnDirs;