%% Kyobi Skutt-Kakaria
% 02.04.2018 - Created
% 02.04.2018 - updated

% this function generates kdes for low to high temp shifts for each of the 5 primary metric
% differences for the paper

function out = generateDiffPlot(data1,error1)


clf
YNames = {'Activity','Turn Bias', 'Turn Direction Switchiness','Wall Proximity','Turn Timing Clumpiness'};
XNames = {'delta z-score','delta z-score', 'delta z-score',...
    'delta z-score','delta z-score'};

for ii = 1:5
    figure
%subplot(2,5,ii)
hold on

idx = any([~isnan(data1(:,ii+10)) ~isnan(data1(:,ii+10))],2);

data1 = data1(idx,:);
error1 = error1(idx,:);


data = data1(:,ii+10) - data1(:,ii+15);
data = (data - nanmean(data))./nanstd(data);
error = sqrt(error1(:,ii+10).^2+error1(:,ii+15).^2)./nanstd(data);

temp = [data1(:,ii+10) data1(:,ii+15)];
dataA = nanmean((temp-nanmean(temp))./nanstd(temp),2);
errorA = nanmean(([error1(:,ii+10), error1(:,ii+15)])./nanstd(temp),2);


dataNull = cell(0);
for kk = 1:100
dataNull{kk} = normrnd(dataA,errorA) - normrnd(dataA,errorA);
end
dataNull = cat(1,dataNull{:});

IBE = bootstrp(1000,@(x,y) IBEV2(x,y),data1(:,[ii+10 ii+15]),...
                error1(:,[ii+10 ii+15]));

            

v1 = nanmean(data)-5*nanstd(data);
v2 = nanmean(data)+5*nanstd(data);

int = ((v2-v1)/500);
bins = (v1+int):int:v2;

hold on
dist1 = fitdist(data,'kernel');
dist2 = fitdist(dataNull,'kernel');

y1 = pdf(dist1,bins);
y2 = pdf(dist2,bins);

hArea1 = area(bins,y2,'LineWidth',2,'FaceAlpha',0.5,'FaceColor',[0.2 0.6 0.2]);
hArea2 = area(bins,y1,'LineWidth',2,'FaceAlpha',0.5,'FaceColor',[0.8 0.2 0.3]);


yMax = max(get(gca,'YLim'));
xMax = max(get(gca,'XLim'));

[~,p1] = ttest2(data,dataNull);
[~,p2] = ztest(nanmean(IBE),0,nanstd(IBE));

vbes = [nanmean(IBE)];
vbee = [nanmean(IBE)];

tStr = strcat({'med_{delta} = '}, num2str(round(nanmedian(data),2)),...
    {'\newlineIBE_{delta} = '}, num2str(round(vbes(1),3)));
hText = text(xMax*0.3, yMax*0.6,tStr,'Interpreter','tex','FontWeight','bold','FontSize',12);

set(gca,'FontSize',16,'FontWeight','bold','XLim',[v1 v2])

if p1 < 0.001
    sym1 = '***';
elseif p1 < 0.01
    sym1 = '**';
elseif p1 < 0.05
    sym1 = '*';
else
    sym1 = '';
end

if p2 < 0.001
    sym2 = '***';
elseif p2 < 0.01
    sym2 = '**';
elseif p2 < 0.05
    sym2 = '*';
else
    sym2 = '';
end

pos1 = hText.Position(1)*2.6;
pos2 = hText.Position(2)*1.03;

text(pos1,pos2,sym1,'FontSize',14,...
    'HorizontalAlignment','center','FontWeight','bold')

text(pos1,pos2*0.94,sym2,'FontSize',14,...
    'HorizontalAlignment','center','FontWeight','bold')

pbaspect([2 1 1])

hLegend = legend([hArea1,hArea2],{'Expected','Observed'},'Location','northwest');

legend('boxoff')

hTitle = title(YNames(ii));

hYLabel = ylabel('Probability Density');
hXLabel = xlabel(XNames(ii));
drawnow

out(ii,:) = [nanmean(IBE) nanstd(IBE) p2];
end
shg
