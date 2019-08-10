centApp = x;

% [a,b] = min(sqrt((centApp(:,1)-shape(:,1)').^2 + (centApp(:,2)-shape(:,2)').^2),[],2);
% wallDistances = [centApp(:,1) - shape(b,1), centApp(:,2) - shape(b,2)];
s = sign(centApp - meanShape);
arm((s(:,1) == -1) & (s(:,2) == 1)) = 1;
arm((s(:,1) == 1) & (s(:,2) == 1)) = 2;
arm(s(:,2) == -1) = 3;

%% Fit each arm independently

clf
ll = 1;
scope = nan(500,3);
lineVals = scope;
params.minDist = 6.5;
params.maxDist = 12;

for ll = 1:3
    idx = (arm == ll) & (distFromCent > params.minDist) & (distFromCent < params.maxDist);
    if ll <= 2
        armFit(ll,:) = polyfit(centApp(idx,1),centApp(idx,2),1);
        scope(:,ll) =  linspace(0,max(centApp(:,1)),500);
        lineVals(:,ll) = polyval(armFit(ll,:) ,scope(:,ll));

    else
        armFit(ll,:) = polyfit(centApp(idx,2),centApp(idx,1),1);
        lineVals(:,ll) =  linspace(0,max(centApp(:,2)),500);
        scope(:,ll) = polyval(armFit(ll,:) ,lineVals(:,ll));
    end
    scatter(centApp(idx,1),centApp(idx,2))
    hold on
end


  


plot(scope,lineVals,'LineWidth',4)
hold on

%% Calculate center more accurately
 m = armFit(:,1);
 b = armFit(:,2);
 
 x = (b(2) - b(1))/(m(1) - m(2));
 y = m(1)*x+b(1);
 
 center = [x,y];
 
s = sign(centApp - center);
arm((s(:,1) == -1) & (s(:,2) == 1)) = 1;
arm((s(:,1) == 1) & (s(:,2) == 1)) = 2;
arm(s(:,2) == -1) = 3;

for ll = 1:3
scatter(centApp(arm==ll,1),centApp(arm==ll,2))
end
%% distances = nan(length(centApp),1);
for ll = 1:3
    [a,b] = min(sqrt((centApp(:,1) - scope(:,ll)').^2 +...
       (centApp(:,2) - lineVals(:,ll)').^2),[],2);
end


%%
clf
distances = nan(length(centApp),1);
for ll = 1:3
    subDat = centApp(arm==ll,:);
    subDist = distFromCent(arm==ll,:);
    [~,b1] = min(subDist); [~,b2] = max(subDist);
    
    if ll <= 2
        armFit = polyfit(subDat(:,1),subDat(:,2),1);
        scope =  linspace(subDat(b2,1),meanShape(1),500);
        lineVals = polyval(armFit,scope);
        
        [a,b] = min(sqrt((subDat(:,1) - scope).^2 +...
            (subDat(:,2) - lineVals).^2),[],2);
    else
        armFit = polyfit(subDat(:,2),subDat(:,1),1);
        scope =  linspace(subDat(b2,2),meanShape(2),500);
        lineVals = polyval(armFit,scope);
        
        [a,b] = min(sqrt((subDat(:,2) - scope).^2 +...
            (subDat(:,1) - lineVals).^2),[],2);
    end
 
    
    if ll == 1
        signs = sign(subDat(:,1) - scope(b)') * -1;
        plot(scope,lineVals,'LineWidth',4)
        hold on
    elseif ll == 2
        signs = sign(subDat(:,1) - scope(b)');
        plot(scope,lineVals,'LineWidth',4)
    else
        signs = sign(subDat(:,1) - lineVals(b)');
        plot(lineVals,scope,'LineWidth',4)
    end
    
    distances(arm==ll) = a.*signs;
    
end



%%
times = turnDirs.tTime{ii};
allTimes = flyTracks.tStamps;
for ll = 1:(length(times)-1)
    idx = find(flyTracks.tStamps > times(ll) & flyTracks.tStamps < times(ll+1));
    m = deltaDist(idx);
    n = distFromCent(idx);
    closePts = cumsum(abs(m),'reverse');
    if max(closePts) > 5
        idx2 = idx(closePts < 15 & closePts > 3);
    else
        idx2 = idx;
    end
    idx3 = idx2(deltaDist(idx2) < -0.5);
    centIdx{ll} = idx3;
    turnIdx{ll} = repmat(turnDirs.tSequence{ii}(ll+1),length(idx3),1);
end
%%
centArr = cat(1,centIdx{:});
turnArr = cat(1,turnIdx{:});
colormap([0.8 0 0.8; 0 0.8 0]);
colorIdx = double(turnArr);
scatter(centApp(centArr,1),centApp(centArr,2),'CData',colorIdx,'Marker','.','MarkerEdgeAlpha',0.5)
drawnow
shg

%% Make fig for fraction R cent versus L cent
dcent = distFromCent(centArr & turnArr);
dwall =  distances(centArr & turnArr);
dcent2 = dcent(dcent>3 & abs(dwall) < 3);
dwall2 = dwall(dcent>3 & abs(dwall) < 3);
%scatter(dcent2,dwall2)



binEdges = 3:1:15;
bins = discretize(dcent2,binEdges);
outVect = nan(length(binEdges),2);
outVectCV = outVect;
for ll = 1:length(binEdges)
    outVect(ll,1) = nanmean(sign(dwall2(bins == ll)));
  %  outVectCV(ll,1) = nanstd(dwall2(bins == ll))/length(dwall2(bins == ll));
end



dcent = distFromCent(centArr & ~turnArr);
dwall =  distances(centArr & ~turnArr);
dcent2 = dcent(dcent>3 & abs(dwall) < 3);
dwall2 = dwall(dcent>3 & abs(dwall) < 3);
%scatter(dcent2,dwall2,'MarkerEdgeColor','red')


binEdges = 3:1:15;
bins = discretize(dcent2,binEdges);
for ll = 1:length(binEdges)
    outVect(ll,2) = nanmean(sign(dwall2(bins == ll)));
  %  outVectCV(ll,2) = nanstd(dwall2(bins == ll))/length(dwall2(bins == ll));
end






%%
centArr = cat(1,centIdx{:});
turnArr = cat(1,turnIdx{:});
colormap([0.8 0 0.8; 0 0.8 0]); 
colorIdx = double(turnArr);
scatter(distFromCent(centArr),distances(centArr),'CData',colorIdx,'Marker','.','MarkerEdgeAlpha',0.5)
drawnow
shg






