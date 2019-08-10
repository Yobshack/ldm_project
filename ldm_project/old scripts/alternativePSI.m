split = 'RH1';
idx = ismember(subArray(:,3),split);
eff = ismember(effectors,'SHI');
effC = ismember(effectors,'ISO');
% 
n = ismember(shorterNames,'PF-LCre');

[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(subArray(:,3)),ngn);

idxC = idx & effC;
idx = idx & eff;

% error1 = subError(idx,:);
% thresh = 0.05;
% idx2 = error1(:,behave) < thresh;

array = subArray(idx,:);
arrayC = subArray(idxC,:);
% array = array(idx2,:);

data1 = subData(idx,:);
data1C = subData(idxC,:);
% data1 = data1(idx2,:);


behave = 2;
% subFlies = randi(length(dataHigh),1000,1);
subFlies = 1:length(data1);
dataHigh = [data1(subFlies,behave+10), data1(subFlies,behave+15)];
dataLow = [data1(subFlies,behave), data1(subFlies,behave+5)];

behave = 2;
% subFlies = randi(length(dataHigh),1000,1);
subFlies = 1:length(data1C);
dataHighC = [data1C(subFlies,behave+10), data1C(subFlies,behave+15)];
dataLowC = [data1C(subFlies,behave), data1C(subFlies,behave+5)];


clf
hAx = axes('OuterPosition',[0.15 0.15 0.5 0.5]);
hold on
hPlot1 = scatter(hAx,dataLow(:,2),dataLow(:,1),'Marker','o',...
    'MarkerEdgeColor','black','MarkerFaceColor',[0 0 0.6],'MarkerFaceAlpha',0.5);

hPlot2 = scatter(hAx,dataHigh(:,2),dataHigh(:,1),'Marker','o',...
    'MarkerEdgeColor','black','MarkerFaceColor',[0.6 0 0],'MarkerFaceAlpha',0.5);

set(gca,'XLim',[0 1],'YLim',[0 1])

hAx2 = axes('OuterPosition',[0.15 0.6 0.5 0.2]);
hold on
edges = 0:0.05:1;
histogram(hAx2,hPlot1.XData,edges,'Normalization','probability',...
    'FaceColor',hPlot1.MarkerFaceColor);
histogram(hAx2,hPlot2.XData,edges,'Normalization','probability',...
    'FaceColor',hPlot2.MarkerFaceColor);

hAx3 = axes('OuterPosition',[0.6 0.15 0.2 0.5]);
view(hAx3,90,270)
hold on
edges = 0:0.05:1;
histogram(hAx3,hPlot1.YData,edges,'Normalization','probability',...
    'FaceColor',hPlot1.MarkerFaceColor);
histogram(hAx3,hPlot2.YData,edges,'Normalization','probability',...
    'FaceColor',hPlot2.MarkerFaceColor);
set([hAx2,hAx3],'XTick','')




  x = hPlot2.XData - hPlot1.XData;
  y = hPlot2.YData - hPlot1.YData;
  
  [theta,rho] = cart2pol(x,y);
  
%   [a,b] = sort(rho,'descend')
%   b = b(~isnan(a));
  [a,b] = sort(abs(dataLow(:,2) - dataLow(:,1)),'descend');
  b = b(~isnan(a));
  a = a(~isnan(a));
  
  top = 1:round(length(b)*0.25);
  b = b(top);
  theta = theta(b);
  rho = rho(b);
  
% for hh = 1:length(b)
%     
%     kk = b(hh);
%     addV = hPlot1.Parent.Position;
%     try
%         hAnno = annotation('arrow',[hPlot1.XData(kk)*addV(3)+addV(1) hPlot2.XData(kk)*addV(3)+addV(1)],...
%             [hPlot1.YData(kk)*addV(4)+addV(2)  hPlot2.YData(kk)*addV(4)+addV(2) ],...
%             'LineWidth',3,'Color',[1 0 0 0.1],'HeadLength',20,'HeadStyle','vback1','HeadWidth',20);
%     end
% end


  x = dataHighC(:,2) - dataLowC(:,2);
  y = dataHighC(:,1) - dataLowC(:,1);
  
  [thetaC,rhoC] = cart2pol(x,y);
  
%   [a,b] = sort(rho,'descend')
%   b = b(~isnan(a));
  [a,b] = sort(abs(dataLowC(:,2) - dataLowC(:,1)),'descend');
  b = b(~isnan(a));
  a = a(~isnan(a));
  
  top = 1:round(length(b)*0.25);
  b = b(top);
  thetaC = thetaC(b);
  rhoC = rhoC(b);

  
  
polaraxes('Position',[0.75 0.4 0.125 0.14])
 

edges = deg2rad(-15:30:345);
  hPol2 = polarhistogram(thetaC,edges,'Normalization','probability');
  hold on
  cntValues = hPol2.Values;
  hPol1 = polarhistogram(theta,edges,'Normalization','probability');
  %hPol1.Values = hPol1.Values - cntValues
   
  
 set(gca,'RTick',[],'ThetaTick',[0:45:360],'ThetaGrid','off')
  
  polaraxes('Position',[0.75 0.15 0.125 0.14])
hold on
  for ii = 1:length(theta)
      hPol = polarplot([0 theta(ii)],[0 rho(ii)],'LineStyle','-','LineWidth',2,'Color',[0 0 0 0.2]);
      hPol = polarplot([0 thetaC(ii)],[0 rhoC(ii)],'LineStyle','-','LineWidth',2,'Color',[0 0.4 0 0.2]);
  end
  
  set(gca,'RTick',[0.25 0.5],'ThetaTick',[0:45:360],'ThetaGrid','off')
  
%   ind1 = theta>0;
%   theta = theta(ind1);
%   rho = rho(ind1);
%   lineAvg = (nansum(abs(theta).*rho))/nansum(rho);
%   
%   polarplot([0 lineAvg], [0 0.5],'LineStyle','-','LineWidth',2,'Color','red')
%   polarplot([0 -(pi - lineAvg)], [0 0.5],'LineStyle','-','LineWidth',2,'Color','red')
%   
outData = [dataLow,dataHigh];
save('LDM_figure2d_1.mat','outData')

%% Individual example

nn = b(5);
time = array{nn,2};
lights = array{nn,11};
turns = array{nn,1};

timeBins = (0:30:180)*60;

tb = nan((length(timeBins)-1)*2,2);
for ii = 1:length(timeBins)-1
    idx = time > timeBins(ii) & time < timeBins(ii+1) & lights;
    tb((ii*2)-1,:) = turnbias(turns(idx));
    idx = time > timeBins(ii) & time < timeBins(ii+1) & ~lights;
    tb(ii*2,:) = turnbias(turns(idx));
end

outData = [tb(1:4,1),tb(9:12,1)];
plot(outData)
line([1 4],[nanmean(tb([1:2:4],1)),nanmean(tb([1:2:4],1))],'color','black')
line([1 4],[nanmean(tb([2:2:4],1)),nanmean(tb([2:2:4],1))],'color','black')

legend('low temp','high temp')

save('LDM_figure2d_2.mat','outData')

%% 

eq1 = @(x) x;

cx = dataHigh(:,2); % dark bias
cy = dataHigh(:,1); % light bias
for ii = 1:size(cx,1)
d = point_to_line([cx(ii),cy(ii),0],[0,0,0],[1,1,0]);

side = sqrt(d/2);

signs = [sign(cx(ii)),sign(cy(ii))];

coordsLine(ii,:) = [cx(ii)+side*(-signs(1)),...
    cy(ii)+side*(-signs(2))]

end
 