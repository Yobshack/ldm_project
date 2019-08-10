%% Kyobi Skutt-Kakaria
% 02.19.2018 - Created
% 02.19.2018 - updated

% this function generates a transition triggered average of IBE while subtracting the pre-transition
% state bias

function hError = generateTurnTMIPlot2(subArray,lightDat,light,color,bin)


null = nan(size(subArray,1),2);
timeBin = 3;
dEnd = [-30:timeBin:5];
tbs = nan(size(subArray,1),2,length(dEnd)-1);

out = nan(length(dEnd),size(subArray,1));
%vect = randi(size(subArray,1),size(subArray,1),1);
for jj = 1:size(subArray,1)

    accumTime = [0;cumsum(lightDat(:,1))];
    turns = subArray{jj,1};
    
    tStamps = subArray{jj,2};
    
    intensity = nan(length(tStamps),1);
    for ii = 1:length(tStamps)
    intensity(ii) = lightDat(sum(tStamps(ii) > accumTime),3);
    end
    ind = intensity == bin;
    
    tStamps = tStamps(ind);
    turns = turns(ind);
    
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
if  length(tStamps) > 50

    mi = nan(length(dEnd),2);
    for ii = 1:(length(dEnd))
        tS = tStamps - dEnd(ii);
        temp = nan(length(tStamps),1);
        for kk = 1:length(tStamps)
            if tS(kk) <= 0
                tS(kk) = 0;
            end
            temp(kk) = lightDat(sum(tS(kk) >= accumTime),2);
            
        end
        if sum(temp) >= 1
        mi(ii,:) = mutualinformation([turns, temp]);
        end
    end
out(:,jj) = mi(:,1);
end
jj
end


alpha = 0.5;
x = dEnd;
y = nanmedian(out,2);
yErr = nanstd(out,[],2)/sqrt(size(out,2));
hError = shadedErrorBar(x,y,yErr);
hError.patch.FaceColor = color;
hError.patch.FaceAlpha = alpha;
% 
% mI = find(dEnd == 0);  
% xFit = dEnd(1:mI);
% fits = fit(xFit',y(1:mI),'exp1');
% coFit = coeffvalues(fits);
% text(-25,y(mI),strcat({'b = '},num2str(round(coFit(2),3))));

% ),'LineStyle','-',...
%     'Marker','o','MarkerFaceColor','auto',...
%     'Color',[0.5 0.5 0.5 alpha],'MarkerSize',5,'CapSize',8,'LineWidth',0.5);

drawnow