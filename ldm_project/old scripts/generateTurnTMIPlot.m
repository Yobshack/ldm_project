%% Kyobi Skutt-Kakaria
% 02.19.2018 - Created
% 02.19.2018 - updated

% this function generates a transition triggered average of IBE while subtracting the pre-transition
% state bias

function [hError1 hError2] = generateTurnTMIPlot(subArray,lightDat,light,color,bin)


null = nan(size(subArray,1),2);
timeBin = 0.5;
dEnd = [-30:timeBin:5];
tbs = nan(size(subArray,1),2,length(dEnd)-1);

out = nan(length(dEnd),size(subArray,1));
out2 = out;
%vect = randi(size(subArray,1),size(subArray,1),1);
for jj = 1:size(subArray,1)

    accumTime = [0;cumsum(lightDat(:,1))];
    turns = subArray{jj,1};
    
    tStamps = subArray{jj,2};
    if size(lightDat,2) == 3
    for ii = 1:length(tStamps);
    intensity(ii) = lightDat(sum(tStamps(ii) > accumTime),3)
    end
    ind = intensity == bin;
    
    tStamps = tStamps(ind);
    turns = turns(ind);
    end
    
    
    temp = nan(length(tStamps),1);
    for kk = 1:length(tStamps)
        temp(kk) = lightDat(sum(tStamps(kk) >= accumTime),2);
    end
    
%     if light
%         tStamps = tStamps(temp);
%         turns = turns(temp);
%     else
%         tStamps = tStamps(~temp);
%         turns = turns(~temp);
%     end
%    
tStampsR = randsample(tStamps,length(tStamps),1);
if  length(tStamps) > 50

    mi = nan(length(dEnd),2);
    nullMi = mi;
    for ii = 1:(length(dEnd))
        tS = tStamps - dEnd(ii);
        temp = nan(length(tStamps),1);
        for kk = 1:length(tStamps)
            if tS(kk) <= 0
                tS(kk) = 0;
            end
            temp(kk) = lightDat(sum(tS(kk) >= accumTime),2);
            
            
        end
        
        mi(ii,:) = mutualinformation([turns, temp]);
        
        % the null is made by randomizing the time stamps of each turn. This retains the structure
        % of the lighting but not the association with a turn.
        tSR = tStampsR - dEnd(ii);
        tempR = nan(length(tStampsR),1);
        for kk = 1:length(tStampsR)
            if tSR(kk) <= 0
                tSR(kk) = 0;
            end
            tempR(kk) = lightDat(sum(tSR(kk) >= accumTime),2);
            
            
        end
        
        nullMi(ii,:) = mutualinformation([turns, tempR]);
    end
out(:,jj) = mi(:,1);
out2(:,jj) = nullMi(:,1);
end
jj
end


alpha = 0.5;
x = dEnd;

hold on
% observed
y = nanmedian(out,2);
yErr = nanstd(out,[],2)/sqrt(size(out,2));
hError1 = shadedErrorBar(x,y,yErr);
hError1.patch.FaceColor = color;
hError1.patch.FaceAlpha = alpha;
hError1.patch.EdgeAlpha = 0;

% null model
y2 = nanmedian(out2,2);
yErr2 = nanstd(out2,[],2)/sqrt(size(out2,2));
hError2 = shadedErrorBar(x,y2,yErr2);
hError2.patch.FaceColor = 'red';
hError2.patch.FaceAlpha = alpha;
hError2.mainLine.Color = [0.5 0 0];
hError2.patch.EdgeAlpha = 0;


% mI = find(dEnd == 0);  
% xFit = dEnd(1:mI);
% fits = fit(xFit',y(1:mI),'exp1');
% coFit = coeffvalues(fits);
% text(-25,y(mI),strcat({'b = '},num2str(round(coFit(2),3))));

% ),'LineStyle','-',...
%     'Marker','o','MarkerFaceColor','auto',...
%     'Color',[0.5 0.5 0.5 alpha],'MarkerSize',5,'CapSize',8,'LineWidth',0.5);

drawnow