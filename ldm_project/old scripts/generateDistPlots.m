%% Kyobi Skutt-Kakaria
% 02.02.2018 - Created
% 02.02.2018 - updated

% this function generates histograms for each of the 5 primary metrics for the paper, very pretty
% plots

function [out,outKernels] = generateDistPlots(data,error,AxH1,AxH2,num,outKernels)


YNames = {'Activity','Turn Bias', 'Turn Direction Switchiness','Wall Proximity','Turn Timing Clumpiness'};
XNames = {'turns per minute','Turn Bias', 'normalized probability changing turn direction',...
    'ipsilateral wall proximity','inter-turn interval variability'};

for ii = 2
% subplot(6,12,num,AxH1)
% %subplot(2,5,ii)


[~,pval1] = ttest(data(:,ii),data(:,ii+5));


idx = ~isnan(data(:,ii)) & ~isnan(data(:,ii+5)) &...
    ~isinf(data(:,ii)) & ~isinf(data(:,ii+5));

v1 = nanmean(data(idx,ii+5))-5*nanstd(data(idx,ii+5));
v2 = nanmean(data(idx,ii+5))+5*nanstd(data(idx,ii+5));

int = ((v2-v1)/500);
bins = 0:0.001:1;

hold on
dist1 = fitdist(data(idx,ii),'kernel');
dist2 = fitdist(data(idx,ii+5),'kernel');

y1 = pdf(dist1,bins);
y2 = pdf(dist2,bins);

dist3 = fitdist(normrnd(mean(data(idx,ii)),error(idx,ii)),'kernel');
dist4 = fitdist(normrnd(mean(data(idx,ii+5)),error(idx,ii+5)),'kernel');

y3 = pdf(dist3,bins);
y4 = pdf(dist4,bins);

outKernels = cat(3,y3,y4);

c1 = [0.2 0.1 0.6];
c2 = [0.8 0.9 0.8];

hArea1 = area(bins,y2,'LineWidth',2,'FaceAlpha',0.75,'FaceColor',c1);
hArea2 = area(bins,y1,'LineWidth',2,'FaceAlpha',0.75,'FaceColor',c2);
% 
% hArea2 = area(bins,y4,'LineWidth',2,'FaceAlpha',0.1,'FaceColor',c1);
% hArea4 = area(bins,y3,'LineWidth',2,'FaceAlpha',0.1,'FaceColor',c2);


yMax = max(get(gca,'YLim'));
xMax = max(get(gca,'XLim'));


rs1 = bootstrp(1000,@(x,y) VBE(x,y),data(:,ii),error(:,ii));
rs2 = bootstrp(1000,@(x,y) VBE(x,y),data(:,ii+5),error(:,ii+5));

[~,pval2] = ttest(rs1,rs2);

vbes = [nanmean(rs1) nanmean(rs2)];
vbee = [nanstd(rs1) nanstd(rs2)];

out1(ii) = (nanmedian(data(:,ii)) - nanmedian(data(:,ii+5)))./nanmedian(data(:,ii+5));
out2(ii) = (vbes(2) - vbes(1))./vbes(2);

out3(ii) = nanmedian(data(:,ii));
out4(ii) = nanmedian(data(:,ii+5));
out5(ii) = vbes(2);
out6(ii) = vbes(1);


% 
% 
% 
% tStr = strcat({'med_{light} = '}, num2str(round(nanmedian(data(:,ii)),2)),...
%     {'\newlinemed_{dark} = '}, num2str(round(nanmedian(data(:,ii+5)),2)),...
%     {'\newlineVBE_{light} = '}, num2str(round(vbes(1),3)),...
%     {'\newlineVBE_{dark} = '}, num2str(round(vbes(2),3)));
% % hText = text(xMax*0.6, yMax*0.6,tStr,'Interpreter','tex','FontWeight','bold','FontSize',12,...
% %     'FontName','Times');

set(gca,'FontSize',24,'FontWeight','bold','XLim',[0 1],'XTick',[0:0.25:1],'TickDir','out')
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontName','Arial','YLim',[0 5],'color','none')
set(gcf,'Renderer','painters')

% if pval1 < 0.001
%     sym1 = '***';
% elseif pval1 < 0.01
%     sym1 = '**';
% elseif pval1 < 0.05
%     sym1 = '*';
% else
%     sym1 = '';
% end
% 
% if pval2 < 0.001
%     sym2 = '***';
% elseif pval2 < 0.01
%     sym2 = '**';
% elseif pval2 < 0.05
%     sym2 = '*';
% else
%     sym2 = '';
% end
% 
% p1 = hText.Position(1)*1.38;
% p2 = hText.Position(2)*1.23;
% hLine1 = line([p1 p1*1.02 p1*1.02 p1],[p2,p2, p2*0.88, p2*0.88],'LineWidth',1.5,'Color','black');
% 
% p1 = hText.Position(1)*1.38;
% p2 = hText.Position(2)*0.96;
% hLine2 = line([p1 p1*1.02 p1*1.02 p1],[p2,p2, p2*0.85, p2*0.85],'LineWidth',1.5,'Color','black');
% 
% text(hLine1.XData(2)*1.04,mean(hLine1.YData(2:3)*0.98),sym1,'FontSize',14,...
%     'HorizontalAlignment','center','FontWeight','bold')
% 
% text(hLine1.XData(2)*1.04,mean(hLine1.YData(2:3)*0.75),sym2,'FontSize',14,...
%     'HorizontalAlignment','center','FontWeight','bold')
% 
% 
hLegend = legend([hArea1,hArea2],{'Dark','Light'},'Location','northwest');
legend('boxoff')
% 
% 
hYLabel = ylabel('Probability Density');
hXLabel = xlabel(XNames(ii));

% 
% f = @(x,y) [mad(x(idx,ii)), mad(normrnd(nanmean(x(idx,ii)),y(idx,ii)))...
%     mad(x(idx,ii+5)), mad(normrnd(nanmean(x(idx,ii+5)),y(idx,ii+5)))];
% barDat = bootstrp(1000,f,data,error);
% 
% subplot(6,12,num,AxH2)
% hold on
% bar(nanmean(barDat),'FaceColor',[0.4 0.4 0.4],'BarWidth',0.75)
% 
% errorbar(nanmean(barDat),nanstd(barDat),'LineStyle','none','Color',[0 0 0],'LineWidth',5,...
%     'CapSize',0)
% names = {'Light Observed','Light Expected','Dark Observed','Dark Expected'};
% set(gca,'FontSize',24,'XTickLabel',[],'XTickLabelRotation',45,'XTick',1:4)


end
shg

out = [out1;out2;out3;out4;out5;out6];

