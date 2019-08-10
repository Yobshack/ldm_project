%% Kyobi Skutt-Kakaria
% 04.01.2018 - Created
% 02.02.2018 - updated

% this function generates histograms for each of the 5 primary metrics for the paper, very pretty
% plots

function [out,out2] = generateDistPlot2(data1,error1)


YNames = {'Activity','Turn Bias', 'Turn Direction Switchiness','Wall Proximity','Turn Timing Clumpiness'};
XNames = {'turns per minute','Turn Bias', 'normalized probability changing turn direction',...
    'ipsilateral wall proximity','inter-turn interval variability'};

for ii = 1:5
    

data = data1(:,ii:5:end);
error = error1(:,ii:5:end);
error = error(all(~isnan(data),2),:);
data = data(all(~isnan(data),2),:);



rs = bootstrp(1000,@(x,y) VBE2(x,y),data,error);


vbes = nanmean(rs(:,1:2));
vbee = nanstd(rs(:,1:2));
% 
% out1(ii) = (nanmedian(data(:,ii)) - nanmedian(data(:,ii+5)))./nanmedian(data(:,ii+5));
% out2(ii) = (vbes(2) - vbes(1))./vbes(2);

out1(ii) = nanmedian(data(:,3))./nanmedian(data(:,1));
out2(ii) = nanmedian(data(:,4))./nanmedian(data(:,2));
out3(ii) = vbes(2);
out4(ii) = vbes(1);
out5(ii,:) = nanmedian(data);
out6(ii,:) = nanmean(rs(:,3:end));



% 
% 
% 
% tStr = strcat({'med_{light} = '}, num2str(round(nanmedian(data(:,ii)),2)),...
%     {'\newlinemed_{dark} = '}, num2str(round(nanmedian(data(:,ii+5)),2)),...
%     {'\newlineVBE_{light} = '}, num2str(round(vbes(1),3)),...
%     {'\newlineVBE_{dark} = '}, num2str(round(vbes(2),3)));
% % hText = text(xMax*0.6, yMax*0.6,tStr,'Interpreter','tex','FontWeight','bold','FontSize',12,...
% %     'FontName','Times');
% 
% set(gca,'FontSize',24,'FontWeight','bold','XLim',[0 1],'XTick',[],'TickDir','out')
% set(gca,'PlotBoxAspectRatio',[1.8 1 1],'FontName','Times','YLim',[0 5],'color','none')
% set(gcf,'Renderer','painters')

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
% hLegend = legend([hArea1,hArea2],{'Dark','Light'},'Location','northwest');
% legend('boxoff')
% % 
% % 
% hYLabel = ylabel('Probability Density');
% hXLabel = xlabel(XNames(ii));
% 
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

out = [out1;out2;out3;out4];
out2 = [out5,out6]';

