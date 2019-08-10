%% Kyobi Skutt-Kakaria
% 02.04.2018 - Created
% 02.04.2018 - updated

% this function generates kdes for low to high temp shifts for each of the 5 primary metric
% differences for the paper


function [out,out2] = generateDiffPlot3(data1,error1)
global cutoff

YNames = {'Activity','Turn Bias', 'Turn Direction Switchiness','Wall Proximity','Turn Timing Clumpiness'};
XNames = {'delta rank','delta rank', 'delta rank',...
    'delta rank','delta rank'};
count = 1;
for ii = 2
%     figure
% %subplot(2,5,ii)
% hold on


data = data1(:,ii:5:end);
error = error1(:,ii:5:end);

temp1 = data(:,1) - data(:,2);
temp2 = mad(temp1)*0.5;

idx = any(data1(:,[1 6]) < 0.02,2);

data(idx,:) = nan;
error(idx,:) = nan;

error = error(all(~isnan(data),2),:);
data = data(all(~isnan(data),2),:);



% [zsort, zsorti] = sort(data);
% 
% for kk = 1:length(zsort)
%     fobs(kk) = diff([find(zsorti(:,1)==kk);find(zsorti(:,2)==kk)])/length(zsort);
% end

    
if size(data,1) < 20
    out(count,:) = [nan,nan];
    out2 = nan;
    continue
end

mean = IBEV12(data,error);
IBE = bootstrp(1000,@(x,y) IBEV12(x,y),data,...
                 error);

% 
% bins = -1:0.005:1;
% 
% % hold on
% dist1 = fitdist(x,'kernel');
% dist2 = fitdist(y,'kernel');
% 
% y1 = pdf(dist1,bins);
% y2 = pdf(dist2,bins);

% hArea1 = area(bins,y2,'LineWidth',2,'FaceAlpha',0.6,'FaceColor',[0.2 0.6 0.2]);
% hArea2 = area(bins,y1,'LineWidth',2,'FaceAlpha',0.6,'FaceColor',[0.8 0.2 0.3]);


% yMax = max(get(gca,'YLim'));
% xMax = max(get(gca,'XLim'));
% 
% % [~,p1] = ttest(data,dataNull);
% [~,p2] = ztest(nanmean(IBE),0,nanstd(IBE));
% % 
% vbes = [nanmean(IBE)];
% % vbee = [nanmean(IBE)];
% % 
% tStr = strcat({'\newlineIBE_{delta} = '}, num2str(round(vbes(1),3)),{' +/- '},...
%     num2str(round(nanstd(IBE),2)));
% hText = text(xMax*0.3, yMax*0.6,tStr,'Interpreter','tex','FontWeight','bold','FontSize',12);

% set(gca,'FontSize',16,'FontWeight','bold','XLim',[v1 v2])
% 
% % if p1 < 0.001
% %     sym1 = '***';
% % elseif p1 < 0.01
% %     sym1 = '**';
% % elseif p1 < 0.05
% %     sym1 = '*';
% % else
% %     sym1 = '';
% % end
% % 
% if p2 < 0.001
%     sym2 = '***';
% elseif p2 < 0.01
%     sym2 = '**';
% elseif p2 < 0.05
%     sym2 = '*';
% else
%     sym2 = '';
% end
% 
% pos1 = hText.Position(1)*2.2;
% pos2 = hText.Position(2)*1.03;
% % 
% % text(p1,p2,sym1,'FontSize',14,...
% %     'HorizontalAlignment','center','FontWeight','bold')
% % 
% text(pos1*1.5,pos2*0.94,sym2,'FontSize',14,...
%     'HorizontalAlignment','center','FontWeight','bold')
% % 
% pbaspect([2 1 1])
% 
% hLegend = legend([hArea1,hArea2],{'Expected','Observed'},'Location','northwest');
% 
% legend('boxoff')
% 
% hTitle = title(YNames(ii));
% 
% hYLabel = ylabel('Probability Density');
% hXLabel = xlabel(XNames(ii));
% drawnow
out(count,:) = [mean(1) nanstd(IBE(:,1),[],1)];
count = count+1;
out2 = size(data,1);
end

