%% Kyobi Skutt-Kakaria
% 02.16.2018 - Created
% 02.16.2018 - updated

% this function generates a transition triggered average of IBE while subtracting the pre-transition
% state bias

function [hPlot,ibes,hErrorR] = generateTTAIBEPlot3(subArray,lightDat,behave,color,light)


null = nan(size(subArray,1),2);
timeBin = 60;
dEnd = [-120:60:360];

tbs = nan(size(subArray,1),2,length(dEnd)-1);
tbsR = tbs;

for jj = 1:size(subArray,1)
    [x1,x2,tt] = calculateRelativeTime(subArray{jj,2},lightDat);
    

    accumTime = [0; cumsum(lightDat(1:(end-1),1))];
    tStamps = subArray{jj,2};
    
    if size(lightDat,2) == 3
    intensity = nan(length(tStamps),1);
    
    for ii = 1:length(tStamps)
    intensity(ii) = lightDat(sum(tStamps(ii) > accumTime),3);
    end
    end
    if light
        ind = x2;
    elseif ~light
        ind = ~x2;
    end
    tStamps = tStamps(ind);
    x = x1(ind);
    walls = subArray{jj,10}(ind);
    y = subArray{jj,1}(ind);
    darkBias = [subArray{jj,13}(:,behave,2,1) subArray{jj,14}(:,behave,2,1) ];
    lightBias = [subArray{jj,13}(:,behave,1,1)  subArray{jj,14}(:,behave,1,1) ];
    
    

        nullD(jj,:) = darkBias;
        nullL(jj,:) = lightBias;



    for ii = 1:(length(dEnd)-1)
        sInd = x >= dEnd(ii) & x < dEnd(ii)+timeBin;
        tS = tStamps(sInd);
        turnVect = y(sInd);
        ws = walls(sInd);
        if length(turnVect)>10
            switch behave
                case 1
                    if light
                        div = (timeBin*sum(tt(:,2)));
                    elseif ~light
                        div = (timeBin*sum(~tt(:,2)));
                    end
                    tbs(jj,:,ii) = [length(turnVect)/div 0.01];
                case 2
                    tbs(jj,:,ii) = turnbias(turnVect);
                case 3
                    tbs(jj,:,ii) = switchiness(turnVect,lightDat,tS);
                case 4
                    tbs(jj,:,ii) = wallDist(ws,turnVect);
                case 5
                    tbs(jj,:,ii) = clumpiness(tS,lightDat);
            end
        end
        
        randX = randsample(x,length(x),1);
        sInd = randX >= dEnd(ii) & randX < dEnd(ii)+timeBin;
        tS = tStamps(sInd);
        turnVect = y(sInd);
        ws = walls(sInd);
        if length(turnVect)>10
            switch behave
                case 1
                    if light
                        div = (timeBin*sum(tt(:,2)));
                    elseif ~light
                        div = (timeBin*sum(~tt(:,2)));
                    end
                    tbsR(jj,:,ii) = [length(turnVect)/div 0.01];
                case 2
                    tbsR(jj,:,ii) = turnbias(turnVect);
                case 3
                    tbsR(jj,:,ii) = switchiness(turnVect,lightDat,tS);
                case 4
                    tbsR(jj,:,ii) = wallDist(ws,turnVect);
                case 5
                    tbsR(jj,:,ii) = clumpiness(tS,lightDat);
            end
        end
    end
    
    
    % plot(dEnd,tbs-nanmean(tbs(1:((length(dEnd)-1)/2))))
    % hold on
end

numRe = 1000;
ibes = nan(numRe,length(dEnd)-1);

    for ii = 1:(length(dEnd)-1)
        data = [nullD(:,1),tbs(:,1,ii)];
        error = [nullD(:,2),tbs(:,2,ii)];
        
        idx = all(~isnan([data error]),2);
        data = data(idx,:);
        error = error(idx,:);
        ibes(:,ii) = bootstrp(numRe,@(x,y) IBEV5(x,y),data,error);
    end

    numRe = 1000;
ibesR = nan(numRe,length(dEnd)-1);

%     for ii = 1:(length(dEnd)-1)
%         data = [null(:,1),tbsR(:,1,ii)];
%         error = [null(:,2),tbsR(:,2,ii)];
%         
%         idx = all(~isnan([data error]),2);
%         data = data(idx,:);
%         error = error(idx,:);
%         ibesR(:,ii) = bootstrp(numRe,@(x,y) IBEV5(x,y),data,error);
%     end
%     

xVect = dEnd(1:(end-1))+diff([dEnd(1),dEnd(2)])/2;

y = abs(squeeze(tbs(:,1,:)) - nullD(:,1));
colorV = abs(nullL(:,1)-nullD(:,1));
%%
clf
hold on
colormap(jet)
sdSet = 0.0006;
bins = 0:0.05:0.6;
dataBins = histc(y,bins);
sds = dataBins.*sdSet;

colorV = colorV(all(~isnan(y),2),:);
y = y(all(~isnan(y),2),:);

cBins = 0:0.05:0.25;

caxis([0 0.25])
colormap(brewermap(length(cBins),'Set1'))
colmap = colormap;
samps = randsample(size(y,1),size(y,1));


for ii = 1:length(cBins)
    if ii < 5
     yM = nanmean(y(colorV > cBins(ii) & colorV < cBins(ii+1),:));
%     yM = interp(yM,100)
    hPlot2 = plot(xVect/60,yM,'Color',colmap(ii,:))
    else
     yM = nanmean(y(colorV > cBins(ii),:));
%     yM = interp(yM,100)
    hPlot2 = plot(xVect/60,yM,'Color',colmap(ii,:))
    end
    set(hPlot2,'LineWidth',2)
end


for ii = 1:300
   fly = samps(ii);
    descY = discretize(y(fly,:),bins);
    x = normrnd(xVect/60,diag(sds(descY,:))');
    hPlot = scatter(x,y(fly,:),'filled','Marker','o','SizeData',75,'CData',...
    repmat(colorV(fly),size(y,2),1),'MarkerFaceAlpha',1,'MarkerEdgeColor','black',...
    'MarkerEdgeAlpha',1,'LineWidth',0.001);
end

shg

%%
% hPlot = plot(xVect,nanmean(ibes));
% set(hPlot(1),'Color','black','LineWidth',1)
% hPatch2 = patch([xVect fliplr(xVect)],...
%     [nanmean(ibes)+nanstd(ibes) fliplr(nanmean(ibes)-nanstd(ibes))],colors{num});
% hPatch2.FaceAlpha = 0.5;
% 
% alpha = 0.5;
% hold on
% hError = errorbar(xVect,nanmean(ibes),nanstd(ibes),'LineStyle',':',...
%     'Marker','o','MarkerFaceColor','auto',...
%     'Color',[color(1,:) alpha],'MarkerSize',6,'CapSize',8,'LineWidth',2);
% 
% % 
% % alpha = 0.1;
% hErrorR = errorbar(xVect,nanmean(ibesR),nanstd(ibesR),'LineStyle','--',...
%     'Marker','o','MarkerFaceColor','auto',...
%     'Color',[1 0 0 alpha],'MarkerSize',10,'CapSize',8,'LineWidth',2);
