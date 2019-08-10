%% Kyobi Skutt-Kakaria
% 02.16.2018 - Created
% 02.16.2018 - updated

% this function generates a transition triggered average of IBE while subtracting the pre-transition
% state bias

function [hPlot,ibes,hErrorR] = generateTTAIBEPlot5(subArray,lightDat,behave,color,light)


lightDat((lightDat(:,3) == 5)|(lightDat(:,3) == 2),3) = 3.5;
lightDat((lightDat(:,3) == 9)|(lightDat(:,3) == 23),3) = 16;
lightDat((lightDat(:,3) == 75)|(lightDat(:,3) == 225),3) = 150;

nullD = nan(size(subArray,1),2);
timeBin = 15;
dEnd = [-30:1:60];

tbs = nan(size(subArray,1),2,length(dEnd)-1,length(unique(lightDat(:,3))));
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
    
    [a,b,c] = unique(intensity);
    for kk = 1:length(a)
        
        ind2 = ind & c == kk;
        sum(ind2);
        
        
        tStamps2 = tStamps(ind2);
        x = x1(ind2);
        walls = subArray{jj,10}(ind2);
        y = subArray{jj,1}(ind2);
        darkBias = turnbias(y(x<0));
        
        
       
        nullD(jj,:,kk) = darkBias;

        
        
        
        for ii = 1:(length(dEnd)-1)
            sInd = x >= dEnd(ii) & x < dEnd(ii)+timeBin;
            tS = tStamps2(sInd);
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
                        tbs(jj,:,ii,kk) = [length(turnVect)/div 0.01];
                    case 2
                        tbs(jj,:,ii,kk) = turnbias(turnVect);
                    case 3
                        tbs(jj,:,ii,kk) = switchiness(turnVect,lightDat,tS);
                    case 4
                        tbs(jj,:,ii,kk) = wallDist(ws,turnVect);
                    case 5
                        tbs(jj,:,ii,kk) = clumpiness(tS,lightDat);
                end
            end
            
        end
        
    end
    % plot(dEnd,tbs-nanmean(tbs(1:((length(dEnd)-1)/2))))
    % hold on
end
%%
numRe = 100;
ibes = nan(numRe,length(dEnd)-1,length(a));
for kk = 1:length(a)
for ii = 1:(length(dEnd)-1)
    if sum(~isnan(tbs(:,1,ii,kk))) > 10
    data = [nullD(:,1,kk),tbs(:,1,ii,kk)];
    error = [nullD(:,2,kk),tbs(:,2,ii,kk)];
    
    idx = all(~isnan([data error]),2);
    data = data(idx,:);
    error = error(idx,:);
    ibes(:,ii,kk) = bootstrp(numRe,@(x,y) mad(diff(x,[],2)),data,error);
    end
end
end



xVect = dEnd(1:(end-1))+diff([dEnd(1),dEnd(2)])/2;
%%
clf
y = squeeze(nanmean(ibes));
e = squeeze(nanstd(ibes));
hold on
for ii = 1:3
errorbar(xVect,y(:,ii),e(:,ii))
end

% set(gca,'YLim',[0.05 0.2])

colorV = abs(nullL(:,1)-nullD(:,1));
%%
clf
hold on
colormap(jet)
for ii = 1:size(y,1)
    x = normrnd(xVect,3)/60;
    hPlot = scatter(x,y(ii,:),'filled','Marker','o','SizeData',50,'CData',...
        repmat(colorV(ii),size(y,2),1),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','black',...
        'MarkerEdgeAlpha',0.5);
end




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
