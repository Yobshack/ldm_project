function [out] = calculateWallDist(turnDirs,flyTracks)

clf;
tic;
for ii = 1:size(flyTracks.centroid,3)

    x = flyTracks.centroid(:,:,ii);
    x = x - min(x);
    
    t = flyTracks.tStamps;
    
    accumulatedTime = [0; cumsum(turnDirs.lightDat(1:(end-1),1))];
    
    lightVect = zeros(length(turnDirs.tTime{ii}),1);
    for hh = 1:length(turnDirs.tTime{ii})
        lightVect(hh) = turnDirs.lightDat(sum(turnDirs.tTime{ii}(hh) - accumulatedTime > 0),2);
    end

    
    % this rotates all the mazes between 65 and 120 so that the analysis works on those
    % mazes as well
    if ii > 64
        x(:,2) = max(x(:,2)) - x(:,2);
        x(:,1) = max(x(:,1)) - x(:,1);
    end
    
    % return maximum yValue
    maxes = max(x(1:100:end,:));
    xmax = maxes(1);
    ymax = maxes(2);
    
    % return shape previously calculated
    xshape =  x(~isnan(x(:,1)),:);
    shape = xshape(boundary(xshape,1),:); 
    if isempty(shape)
        continue
    end  
    
    % filter messed up mazes 
    % this could maybe be revised
    fullArea = xmax*ymax;
    shapeArea = polyarea(shape(:,1),shape(:,2));
    proArea(ii) = shapeArea/fullArea;

    
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
    approaches = ((distFromCent < farThresh) & (distFromCent > closeThresh)) &...
        (deltaDist <= -0.05);
    
    % find turn indexes
    tIdx = find(~isnan(flyTracks.rightTurns(:,ii)));
    tIdx = tIdx(2:end);
    
    % the meat of this script the idea is to calculate which wall the animal is closer to in every
    % approach frame
    avgWallDist = nan(length(turnDirs.tSequence{ii}),1);
    avgCentDist = avgWallDist;
    allWallDist = cell(0);
    allCentDist = allWallDist;
    allCentApp = allWallDist;
    allTurnDir = allWallDist;
    allLight = allWallDist;
    allTStamps = allWallDist;
    for jj = 1:length(turnDirs.tSequence{ii})
        
        % finds frames that are approach frames and also fall within a certain amount of time before
        % the fly passes through the choice point
        b1 = tIdx(jj);
        distBack = zeros(1,1);
        count = 1;
        while distBack < farThresh
            if count == (b1 - 1)
              distBack = farThresh;
            else
            distBack = distBack + abs(deltaDist(b1-count));
            count = count+1;
            end
        end
        
        b2 = b1-count;
        cent = x(b2:b1,:);
        centApp = cent(approaches(b2:b1),:);
        
        tT = t(b2:b1,:);
        tApp = tT(approaches(b2:b1),:);
        
        distCent = distFromCent(b2:b1,:);
        distCentApp = distCent(approaches(b2:b1,:));
        
        
        % for each centroid, this finds a point on the boundary outline that is closest to the
        % centroid at that frame and calculates distance, first column is right, second is left
        
        wallDist = nan(size(centApp,1),2);
        for ll = 1:size(centApp,1)
            armlogical = centApp(ll,:) - meanShape > 0;
            shapelogical = shape - meanShape > 0;
            if ~armlogical(1) && armlogical(2)
                s = shape(~shapelogical(:,1) & shapelogical(:,2),:);
                dist = centApp(ll,:) - s;
                l = min(sqrt(sum(dist(dist(:,2) <= 0,:).^2,2)));
                r = min(sqrt(sum(dist(dist(:,2) >= 0,:).^2,2)));               
            elseif armlogical(1) && armlogical(2)
                s = shape(shapelogical(:,1) & shapelogical(:,2),:);
                dist = centApp(ll,:) - s;
                l = min(sqrt(sum(dist(dist(:,2) >= 0,:).^2,2)));
                r = min(sqrt(sum(dist(dist(:,2) <= 0,:).^2,2)));   
            elseif ~armlogical(2)
                s = shape(~shapelogical(:,2),:);
                dist = centApp(ll,:) - s;
                l = min(sqrt(sum(dist(dist(:,1) >= 0,:).^2,2)));
                r = min(sqrt(sum(dist(dist(:,1) <= 0,:).^2,2)));
            else
                l = nan;
                r = nan;
            end
            if isempty(l) || isempty(r)
                wallDist(ll,:) = [nan nan];
            else
                wallDist(ll,:) = [l,r];
            end
        end
        
        % this calculates the proportional position between the two sides of the maze arms
        proportionWD = wallDist(:,2)./sum(wallDist,2);
        avgWallDist(jj) = nanmean(proportionWD);
        avgCentDist(jj) = nanmean(distCentApp);
        
        allWallDist{jj} = proportionWD;
        allCentDist{jj} = distCentApp;
        allCentApp{jj} = centApp;
        allTurnDir{jj} = repmat(turnDirs.tSequence{ii}(jj),size(centApp,1),1);
        allTStamps{jj} = tApp;

        allLight{jj} = repmat(lightVect(jj),size(centApp,1),1);

               
    end
    
    turnDirs.avgWallDist{ii} = avgWallDist;
    turnDirs.avgCentDist{ii} = avgCentDist;
    turnDirs.avgLight{ii} = lightVect;
    turnDirs.allWallDist{ii} = cat(1,allWallDist{:});
    turnDirs.allCentDist{ii} = cat(1,allCentDist{:});
    turnDirs.allCentApp{ii} = cat(1,allCentApp{:});
    turnDirs.allTurnDir{ii} = cat(1,allTurnDir{:});
    turnDirs.allLight{ii} = cat(1,allLight{:});
    turnDirs.allTStamps{ii} = cat(1,allTStamps{:});
   
    
end
toc;
turnDirs.proArea = proArea;
out = turnDirs;

% %% example figs
% 
% clf
% fly = 34
% 
% colors = [1 0 1;0 1 0];
% scatter(turnDirs.allCentApp{fly}(:,1),turnDirs.allCentApp{fly}(:,2),...
%     [],double(turnDirs.allTurnDir{fly}),'MarkerEdgeAlpha',0.4,'Marker','.')
% colormap(gca,colors)
% pbaspect([1 1 1])
% shg
% %%
% clf
% clear tData lIdx rIdx wData lData rData centData centApp
% for kk = 1:2
%     for jj = 1:120
%         
%         fly = jj;
%         if kk == 1
%             lightstatus = logical(turnDirs.allLight{fly});
%         else
%             lightstatus = ~logical(turnDirs.allLight{fly});
%         end
%         
%         tData{jj,kk} = turnDirs.allTurnDir{fly}(lightstatus);
%         lIdx = ~tData{jj,kk};
%         rIdx = tData{jj,kk};
%         wData{jj,kk} = turnDirs.allWallDist{fly}(lightstatus);
%         lData{jj,kk} = wData{jj,kk}(lIdx);
%         rData{jj,kk} = wData{jj,kk}(rIdx);
%         nanmean(rData{jj,kk}) - nanmean(lData{jj,kk});
%         centData{jj,kk} = turnDirs.allCentDist{fly}(lightstatus);
%         centApp{jj,kk} = turnDirs.allCentApp{fly}(lightstatus,:);
%     end
%     ldmWallVect(:,kk) = cellfun(@nanmean,rData(:,kk)) - cellfun(@nanmean,lData(:,kk));
% end
% [a,b] = sort(ldmWallVect(:,2) - ldmWallVect(:,1))
% 
% %%
% fly = 27;
% lightStatus = 2;
% tData1 = logical(cat(1,tData{fly,lightStatus}));
% lIdx = ~tData1;
% rIdx = tData1;
% wData1 = cat(1,wData{fly,lightStatus});
% lData1 = cat(1,lData{fly,lightStatus});
% rData1 = cat(1,rData{fly,lightStatus});
% centData1 = cat(1,centData{fly,lightStatus});
% centApp1 = cat(1,centApp{fly,lightStatus});
% 
% 
% subplot(3,4,[1,2,5,6,9,10])
% 
% scatter(centApp1(:,1),centApp1(:,2),...
%     [],double(tData1),'MarkerEdgeAlpha',0.3,'Marker','.')
% colormap(gca,colors)
% pbaspect([1 1 1])
% 
% binEdges = 4:2:14;
% 
% b1 = discretize(centData1(lIdx),binEdges);
% b2 = discretize(centData1(rIdx),binEdges);
% 
% binVect = nan(length(binEdges),2);
% binCV = binVect;
% binCell = cell(0);
% for ii = 1:max(b1)
%     binVect(ii,1) = nanmean(lData1(b1==ii));
%     binVect(ii,2) = nanmean(rData1(b2==ii));
%     binCV(ii,1) = nanstd(lData1(b1==ii));
%     binCV(ii,2) = nanstd(rData1(b2==ii));
%     binCell{ii}{1} = lData1(b1==ii);
%     binCell{ii}{2} = rData1(b2==ii);
% end
% 
% binVect = [smooth(binVect(:,1),10), smooth(binVect(:,2),10)];
% 
% subplot(3,4,[3 7 11])
% scatter(wData1,centData1,...
%     'CData',double(tData1),'MarkerEdgeAlpha',0.5,'Marker','.')
% colormap(gca,colors)
% hold on
% pbaspect([1,3,1])
% 
% % 
% % plot(-binEdges',binVect(:,1),'-o','LineWidth',3,'Color',[0.3 0 0.3],'Marker','none')
% % plot(-binEdges',binVect(:,2),'-o','LineWidth',3,'Color',[0 0.3 0],'Marker','none')
% % colormap(gca,colors)
% % line(-[3 20],[nanmean(binVect(:)) nanmean(binVect(:))],'color','black')
% 
% 
% histBinEdges = 0:0.05:1;
% 
% subplot(3,4,4)
% 
% cellNum = max(b1);
% h = histogram(binCell{cellNum}{1},histBinEdges,'Normalization','probability','FaceColor',colors(1,:))
% hold on
% line([nanmedian(binCell{cellNum}{1}) nanmedian(binCell{cellNum}{1})], [0 0.25],...
%     'color','black','LineWidth',2)
% histogram(binCell{cellNum}{2},histBinEdges,'Normalization','probability','FaceColor',colors(2,:))
% ylim([0 0.3])
% line([nanmedian(binCell{cellNum}{2}) nanmedian(binCell{cellNum}{2})], [0 0.25],...
%     'color','black','LineWidth',2)
% 
% subplot(3,4,8)
% 
% cellNum = floor(max(b1)/2);
% histogram(binCell{cellNum}{1},histBinEdges,...
%     'Normalization','probability','FaceColor',colors(1,:))
% hold on
% line([nanmedian(binCell{cellNum}{1}) nanmedian(binCell{cellNum}{1})], [0 0.25],...
%     'color','black','LineWidth',2)
% histogram(binCell{cellNum}{2},histBinEdges,'Normalization','probability','FaceColor',colors(2,:))
% ylim([0 0.3])
% line([nanmedian(binCell{cellNum}{2}) nanmedian(binCell{cellNum}{2})], [0 0.25],...
%     'color','black','LineWidth',2)
% subplot(3,4,12)
% 
% 
% cellNum = 1;
% histogram(binCell{cellNum}{1},histBinEdges,'Normalization','probability','FaceColor',colors(1,:))
% hold on
% line([nanmedian(binCell{cellNum}{1}) nanmedian(binCell{cellNum}{1})], [0 0.25],...
%     'color','black','LineWidth',2)
% histogram(binCell{cellNum}{2},histBinEdges,'Normalization','probability','FaceColor',colors(2,:))
% ylim([0 0.3])
% line([nanmedian(binCell{cellNum}{2}) nanmedian(binCell{cellNum}{2})], [0 0.25],...
%     'color','black','LineWidth',2)
% shg
% 
% % % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
