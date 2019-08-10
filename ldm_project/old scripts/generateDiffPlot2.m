%% Kyobi Skutt-Kakaria
% 02.04.2018 - Created
% 02.04.2018 - updated

% this function generates kdes for low to high temp shifts for each of the 5 primary metric
% differences for the paper

function out = generateDiffPlot2(data1,error1)


YNames = {'Activity','Turn Bias', 'Turn Direction Switchiness','Wall Proximity','Turn Timing Clumpiness'};
XNames = {'delta rank','delta rank', 'delta rank',...
    'delta rank','delta rank'};
count = 1;
for ii = 1:5
%     figure
% %subplot(2,5,ii)
% hold on

d1 = data1(:,ii);
d2 = data1(:,ii+5);

e1 = error1(:,ii);
e2 = error1(:,ii+5);


idx = ~isnan(d1) & ~isnan(d2);

data = [d1(idx), d2(idx)];

error = [e1(idx),e2(idx)];

% [zsort, zsorti] = sort(data);
% 
% for kk = 1:length(zsort)
%     fobs(kk) = diff([find(zsorti(:,1)==kk);find(zsorti(:,2)==kk)])/length(zsort);
% end

dataNull = cell(0);
for kk = 1:100
dataNull{kk} = [normrnd(nanmean(data,2),nanmean(error,2)),normrnd(nanmean(data,2),nanmean(error,2))];
end
dataNull = cat(1,dataNull{:});

% [sv,si] = sort(data);
% [sve,svi] = sort(dataNull);
% 
% x = arrayfun(@(x) find(si(:,1) == x) - find(si(:,2) == x),1:length(si))/length(si);
% 
% y = arrayfun(@(y) find(svi(:,1) == y) - find(svi(:,2) == y),1:length(svi))/length(svi);

    x = data(:,1)-data(:,2);
    y = dataNull(:,1) - dataNull(:,2);
    
if isempty(data)
continue
end

IBE = bootstrp(1000,@(x,y) IBEV3(x,y),data,...
                 error);


bins = -1:0.005:1;

% hold on
dist1 = fitdist(x,'kernel');
dist2 = fitdist(y,'kernel');

y1 = pdf(dist1,bins);
y2 = pdf(dist2,bins);

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
out(count,:) = [nanmean(IBE) nanstd(IBE)];
count = count+1;
end

