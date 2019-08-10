%% Kyobi Skutt-Kakaria
% 02.16.2018 - Created
% 02.16.2018 - updated

% this function generates a transition triggered average of IBE while subtracting the pre-transition
% state bias

function hError = generateTTAPlot3(subArray,lightDat,behave,color,light)


null = nan(size(subArray,1),2);
timeBin = 60;
dEnd = [-400:30:400];
tbs = nan(size(subArray,1),2,length(dEnd)-1);

for jj = 1:size(subArray,1)
    [x1,x2,tt] = calculateRelativeTime(subArray{jj,2},lightDat);
    

    accumTime = [0; cumsum(lightDat(1:(end-1),1))];
    tStamps = subArray{jj,2};
    intensity = nan(length(tStamps),1);
    
    for ii = 1:length(tStamps)
    intensity(ii) = lightDat(sum(tStamps(ii) > accumTime),3);
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
    
    
    if light
        null(jj,:) = darkBias;
    elseif ~light
        null(jj,:) = lightBias;
    end

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
    end
    
    
    % plot(dEnd,tbs-nanmean(tbs(1:((length(dEnd)-1)/2))))
    % hold on
end

numRe = 1000;
ibes = nan(numRe,length(dEnd)-1);

    for ii = 1:(length(dEnd)-1)
        data = [null(:,1),tbs(:,1,ii)];
        error = [null(:,2),tbs(:,2,ii)];
        
        idx = all(~isnan([data error]),2);
        data = data(idx,:);
        error = error(idx,:);
        ibes(:,ii) = bootstrp(numRe,@(x,y) IBEV3(x,y),data,error);
    end


xVect = dEnd(1:(end-1))+diff([dEnd(1),dEnd(2)])/2;
% hPlot = plot(xVect,nanmean(ibes));
% set(hPlot(1),'Color','black','LineWidth',1)
% hPatch2 = patch([xVect fliplr(xVect)],...
%     [nanmean(ibes)+nanstd(ibes) fliplr(nanmean(ibes)-nanstd(ibes))],colors{num});
% hPatch2.FaceAlpha = 0.5;

alpha = 0.5;
hError = errorbar(xVect,nanmean(ibes),nanstd(ibes),'LineStyle','none',...
    'Marker','o','MarkerFaceColor','auto',...
    'Color',[color alpha],'MarkerSize',10,'CapSize',8,'LineWidth',2);

drawnow