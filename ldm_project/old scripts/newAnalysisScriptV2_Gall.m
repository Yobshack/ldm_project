%% Select files to analyze

% opens a ui to select mat files to analyze
fullPath = pwd;
files = uigetfile('*.mat','MultiSelect','on');
if iscell(files) ~= 1
    files = {files};
end

% loads an experimental design file, script expects a n x 3 array where n is the number of bins in
% the stimulus sequence. col 1 is the bin length in seconds, col 2 is a logical vector indicating
% whether the `lights were on or off in that bin. col 3 is a vector of pwm values that reflect the
% pwm setting of the arduino controlling the led light board.
lightDat = dlmread('lightSequence_15Min.txt');
lightDat = [lightDat;lightDat(end,:)];
lightDat(end,1) = lightDat(end,1)*10;
tic
%turnDirs = batch(clust,@extractTurns,1,{files,lightDat});
toc

%turnDirs = cell(0);
tic
for ii = 1:length(files)
    load(files{ii}) 
    turnDirs = extractTurns(flyTracks,files{ii},lightDat);
    toc
    save(strcat(files{ii}(1:(end-4)),'_processed.mat'),'turnDirs')
end
toc

%% Make large cell array to handle all data that I can query for groups
% 
% 
% % files = uigetfile('*.mat','MultiSelect','on');
% tmpC = cell(0);
% for ii = 1:length(files)
%     load(files{ii})
%     ii
%     for ll = 1:length(turnDirs.tLabels(:,1))
%         if ~isempty(turnDirs.tLabels{ll,1})
%         turnDirs.tLabels{ll,1} = strcat(turnDirs.tLabels{ll,1},'_Brp');
%         end
%     end
%     % need to heal a few of the labels files because the number of roi's was not 120 on some trays
%     if size(turnDirs.tLabels,1) > size(turnDirs.tSequence,1)
%         turnDirs.tLabels = turnDirs.tLabels(1:size(turnDirs.tSequence,1),:);
%     elseif size(turnDirs.tLabels,1) < size(turnDirs.tSequence,1)
%         turnDirs.tLabels = [turnDirs.tLabels;repmat(turnDirs.tLabels(end,:),...
%             size(turnDirs.tSequence,1)-size(turnDirs.tLabels,1))];
%     end     
%     
%     if istable(turnDirs.tLabels)
%         turnDirs.tLabels = table2cell(turnDirs.tLabels);
%     end
%         
%     
%     
%     dataCell = cell(size(turnDirs.tSequence,1),12);
%     for kk = 1:size(turnDirs.tSequence,1)
%         if length(turnDirs.avgWallDist) < kk || isempty(turnDirs.tLabels{kk,1})
%             continue
%         end
%         % separate geno and effector
%         labs = turnDirs.tLabels(kk,:);
%         labs(6) = cell(1);
%         labsSplit = strsplit(char(labs(1)),'_');
%         labs = [labsSplit labs(2:(end-1))];
%         dataCell(kk,:) = [turnDirs.tSequence(kk),turnDirs.tTime(kk),labs,...
%             {turnDirs.dateAndTime},turnDirs.avgWallDist(kk),...
%             {logical(turnDirs.avgLight{kk})},turnDirs.proArea(kk)];
%     end
%     tmpC{ii} = dataCell;
% end
% 
% cellColNames = {'turnSequence','turnTimes','genotype','effector','sex',...
%     'treatment','box','maze','dateAndTime','wallDistance','lightStatus',...
%     'totalMazeArea'};
% 
% % concatenate all data into the same cell array
% allDataCell = cat(1,tmpC{:});
% 
% % take only non-empty cells
% allDataCell = allDataCell(~sum(cellfun(@isempty,allDataCell),2) > 0,:);
% 
%     % filter for proper maze size
% a = cell2mat(allDataCell(:,12));
% idx = a > (median(a) - 2*mad(a)) & a < (median(a) + 2*mad(a));
% 
% allDataCell = allDataCell(idx,:);

%%
files = uigetfile('*.mat','MultiSelect','on');
[allDataCell, cellColNames] = parseProcessedFiles(files);

%%
importAnno = GetGoogleSpreadsheet('1imo17lojEJHSthqFs9pHBWPjG_F-CRBeTUWztrszt0A');

%importAnno([139, 147, 152, 173],:) = []%importAnno([140, 148, 153, 170, 172],[10 9]) 
%importAnno([107,152],:) = [];
importAnno1 = importAnno(101:end,:);
geno = 'SS02252';
idx = strcmp(allDataCell(:,3),geno);
idx = logical(ones(size(allDataCell,1),1));
adca = allDataCell(idx,:);

importAnno1 = importAnno;
adca = allDataCell;

arrayDayMonth = nan(size(adca,1),4);
for ii = 1:size(adca,1)
    if ~isdatetime(adca{ii,9})
    adca{ii,9} = datetime(strcat(adca{ii,9}(3),adca{ii,9}(1),adca{ii,9}(2),...
    adca{ii,9}(4),adca{ii,9}(5)),'InputFormat',...
    'yyyyMMddHHmm');
    end
    arrayDayMonth(ii,:) = [month(adca{ii,9}) day(adca{ii,9}) adca{ii,8} adca{ii,7}];
end


annoDayMonth = nan(size(importAnno1,1),4);
for ii = 1:size(importAnno1,1)
    if ~isdatetime(importAnno1{ii,3})
    importAnno1(ii,3) = {datetime(char(importAnno1{ii,3}),'InputFormat','MM_dd_yyyy')};
    end
    annoDayMonth(ii,:) = [month(importAnno1{ii,3}), day(importAnno1{ii,3}), importAnno1{ii,1},...
        str2double(importAnno1{ii,15})];
end

[a,b,c] = intersect(arrayDayMonth,annoDayMonth,'rows');

adca = adca(b,:);
importAnno1 = importAnno1(c,:);

idx = cellfun(@isempty,importAnno1(:,18)) | cellfun(@isempty,importAnno1(:,19));
adca(idx,:) = [];
importAnno1(idx,:) = [];

%%

turnFilter = 50;

% cell fills matrix as follows, animal x behavior x light x timebin x intensity bin
allData = cell(0); allError = allData;
times = {[0 180]*60};
bins = 120*60;
for hh = 1:size(adca,1)
    hh
    [allData{hh}, allError{hh}] = generateBehavioralMetrics3(adca(hh,:),lightDat,times,...
        bins,turnFilter);
end
    %cellColNames(13:14) = {'LightLow','DarkLow','LightHigh','DarkHigh'};
    
    
array = adca;

array(:,13) = allData;
array(:,14) = allError;

array(:,3) = upper(array(:,3));

data = cat(1,array{:,13});
error = cat(1,array{:,14});



%%


% Change maze labeling for each new maze id
thresh = 7.5;

idx = cell2mat(importAnno1(:,18)) > thresh & cell2mat(importAnno1(:,19)) > thresh ...
    & ~any(cellfun(@isempty,importAnno1(:,[16 17])),2) & all(data(:,[1 6]) > 0.01,2);
%idx =  ~any(cellfun(@isempty,importAnno1(:,[9 10])),2) & all(data(:,[1 6]) > 0.02,2);

importAnno2 = importAnno1(idx,:);
array1 = array(idx,:,:);
data1 = data(idx,:,:);
error1 = error(idx,:,:);

% f = @(x) (x-0.5)./nanstd(x,[],1);
% %f = @(x) x*1
% 
% data1 = f(data1);

data2 = [data1(:,:,1) data1(:,:,2) data1(:,:,1)-data1(:,:,2)];


opDat = cell2mat((importAnno2(:,16:17)));
% 
% flipV = [9 13 16 24 26]
% opDat(flipV,:) = opDat(flipV,[2,1]);

corrMat1 = [sum([opDat(:,1),opDat(:,2)],2) log2(opDat(:,1)./opDat(:,2)) data2];

figure
scatter(corrMat1(:,14),corrMat1(:,1),'filled');

%text(corrMat1(:,1),corrMat1(:,14),num2str(a))
corr(corrMat1(:,14),corrMat1(:,1),'rows','pairwise','type','Spearman')

figure
scatter(corrMat1(:,14),corrMat1(:,2),'filled');

name = strcat(importAnno2(:,4),'_',importAnno2(:,1));
text(corrMat1(:,14),corrMat1(:,2),name,'Interpreter','none')
reSamp = bootstrp(10000,@(x) corr(x(:,2),x(:,14),'rows','pairwise','type','Spearman'),corrMat1);
mean(reSamp)
shg


%%
[~,p] = corrplot(corrMat1(:,[ 1 2 4 9 14 6 11 16]),'testR','on','type','Spearman')


%%
figure
bar(nanmean(opDat),'BarWidth',0.75,'FaceColor',[0.4 0.4 0.4])
hold on
reSamp = bootstrp(1000,@nanmean,opDat)
errorbar(nanmean(opDat),nanstd(reSamp),'Color','black','LineWidth',5,'LineStyle','none','CapSize',0)
set(gca,'XLim',[0 3],'XTickLabel',{'L','R'},'FontSize',36,'FontWeight','bold','FontName','times')
ylabel('Total Bouton Volume (um3)')
pbaspect([1 1.5 1])
%set(gcf,'PaperSize',[1 1])
set(gca,'YLim',[0 8000])
print(strcat(geno,'_Total'),'-dpdf','-fillpage')


%%
figure

hold on
for ii = 1:size(opDat,1)
    plot([1 2],opDat(ii,:),'Marker','o','MarkerFaceColor','auto','Color',[0 0 0],'MarkerSize',15)
end

set(gca,'XLim',[0.5 2.5],'XTickLabel',{'L','R'},'FontSize',36,'XTick',[1:2],...
    'FontWeight','bold','FontName','times')
pbaspect([1 1.5 1])
ylabel('Total Bouton Volume (um3)')
%set(gcf,'PaperSize',[1 1])
set(gca,'YLim',[0 1400 ])
print(strcat(geno,'_Pairs'),'-dpdf','-fillpage')
shg

%% 
figure


scatter(corrMat1(:,14),corrMat1(:,2),'filled','SizeData',repmat(200,size(corrMat1,1),1),...
    'MarkerEdgeColor','black');

set(gca,'XLim',[-0.5 0.5],'YLim',[-1.5 1.5],'FontSize',32,'FontWeight','bold','FontName','arial')
ylabel('Log2(LVol/RVol)')
xlabel('TBlight - TBDark')
pbaspect([1 1 1])
[corr1,p] = corrcoef(corrMat1(:,14),corrMat1(:,2),'rows','pairwise')
%  corr1 = corr(corrMat1(:,14),corrMat1(:,2),'rows','pairwise','type','Spearman')
%  p = 1;
hLS = lsline
set(hLS,'LineWidth',1,'Color',[0 0 0])
hText1 = text(0.05,-0.4,strcat({'r = '},num2str(round(corr1(1,2),2))));
hText2 = text(0.05,-0.6,[strcat({'p = '},num2str(round(p(1,2),2))),strcat({'n = '},...
    num2str(size(corrMat1,1)))] );
set([hText1,hText2],'FontSize',24)
set(gcf,'Renderer','painters')
%set(gcf,'PaperSize',[1 1])
print(strcat(geno,'_Corr'),'-dpdf','-fillpage')


shg


%% 
figure
scatter(corrMat1(:,9),corrMat1(:,2),'filled','SizeData',repmat(200,size(corrMat1,1),1),...
    'MarkerEdgeColor','black');

set(gca,'XLim',[0 1],'YLim',[-1.5 1.5],'FontSize',32,'FontWeight','bold','FontName','arial')
ylabel('Log2(LVol/RVol)')
xlabel('TBDark')
pbaspect([1 1 1])
[corr1,p] = corrcoef(corrMat1(:,9),corrMat1(:,2),'rows','pairwise')
%  corr1 = corr(corrMat1(:,14),corrMat1(:,2),'rows','pairwise','type','Spearman')
%  p = 1;
hLS = lsline
set(hLS,'LineWidth',1,'Color',[0 0 0])
hText1 = text(0.05,-0.4,strcat({'r = '},num2str(round(corr1(1,2),3))));
hText2 = text(0.05,-0.6,strcat({'p = '},num2str(round(p(1,2),2))));
set([hText1,hText2],'FontSize',24)
set(gcf,'Renderer','painters')
%set(gcf,'PaperSize',[1 1])
print(strcat(geno,'_Dark_Corr'),'-dpdf','-fillpage')


shg
%%
figure
scatter(corrMat1(:,4),corrMat1(:,2),'filled','SizeData',repmat(200,size(corrMat1,1),1),...
    'MarkerEdgeColor','black');

set(gca,'XLim',[0 1],'YLim',[-1.5 1.5],'FontSize',32,'FontWeight','bold','FontName','arial')
ylabel('Log2(LVol/RVol)')
xlabel('TBlight')
pbaspect([1 1 1])
[corr1,p] = corrcoef(corrMat1(:,4),corrMat1(:,2),'rows','pairwise')
%  corr1 = corr(corrMat1(:,14),corrMat1(:,2),'rows','pairwise','type','Spearman')
%  p = 1;
hLS = lsline
set(hLS,'LineWidth',1,'Color',[0 0 0])
hText1 = text(0.05,-0.4,strcat({'r = '},num2str(round(corr1(1,2),2))));
hText2 = text(0.05,-0.6,strcat({'p = '},num2str(round(p(1,2),2))));
set([hText1,hText2],'FontSize',24)
set(gcf,'Renderer','painters')
%set(gcf,'PaperSize',[1 1])
print(strcat(geno,'_Light_Corr'),'-dpdf','-fillpage')


shg

%%
figure
scatter(corrMat1(:,14),corrMat1(:,1),'filled','SizeData',repmat(200,size(corrMat1,1),1),...
    'MarkerEdgeColor','black');

set(gca,'XLim',[-0.25 0.25],'FontSize',32,'FontWeight','bold','FontName','Arial','YLim',[0 2000])
ylabel('LVol + RVol')
xlabel('TBlight - TBDark')
pbaspect([1 1 1])
[corr1,p] = corrcoef(corrMat1(:,14),corrMat1(:,1),'rows','pairwise')
%  corr1 = corr(corrMat1(:,14),corrMat1(:,2),'rows','pairwise','type','Spearman')
%  p = 1;
hLS = lsline
set(hLS,'LineWidth',1,'Color',[0 0 0])
hText1 = text(0.05,2000,strcat({'r = '},num2str(round(corr1(1,2),2))));
hText2 = text(0.05,800,strcat({'p = '},num2str(round(p(1,2),2))));
set([hText1,hText2],'FontSize',24)
set(gcf,'Renderer','painters')
%set(gcf,'PaperSize',[1 1])
print(strcat(geno,'_Vol_Corr'),'-dpdf','-fillpage')


shg

%% Generate a transition triggered average

figure
colors = {[0 0 0],[0.8 0.2 0],[0 0.6 0.2]}; 
alpha = 0.5;

behave = 2;
light = 1;

c1 = [0.3 0.18 0.7];
c2 = [0.4 0.45 0.5];

subArray = array;

hold on
%yyaxis left
hError1 = generateTTAIBEPlot3(subArray,lightDat,behave,[0 0 0],light);
% figure
% hold on
% set(gca,'FontSize',42,'YColor',[0 0 0],'YLim',[0.03 0.13],'YTick',[0:0.025:0.15],'XTick',...
%     -120:120:360,'FontWeight','bold','FontName','times','XLim',[-120 240],'TickDir','out')

    pbaspect([1.5 1 1])
xlabel('Time (s)')
set(gcf,'Renderer','painters')
%ylabel('Activity Corrected LDM')
%line([0 0],[0 0.125],'Color','black','LineStyle','--','LineWidth',3)
if light
    hPatch = patch([-120 -120 0 0],[0 0.15 0.15 0],[0 0 0]);

else
    hPatch = patch([0 0 360 360],[0 0.15 0.15 0],[0 0 0]) ;
end
set(hPatch,'FaceAlpha',0.15)
%%
clf
generateTTALdmVersusMedWall2(subArray,lightDat,0,corrMat1(:,2))
%% Mean Shift plots


[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);

idx = 1:size(array,1)

[vars1,vars2,vars3,vars4,vars5,vars6] = generateDistPlots(data(idx,:),error(idx,:));


    

%% 
figure
imdata = log2(rdivide(medianMatrixL(:,:,2),medianMatrixL(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Light')
savefig('ShiMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = log2(rdivide(medianMatrixD(:,:,2),medianMatrixD(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Dark')
savefig('ShiMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = log2(rdivide(medianMatrixL(:,:,3),medianMatrixL(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESTrp - ESIso')
title('Trp Mean Light')
savefig('ShiMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = log2(rdivide(medianMatrixD(:,:,3),medianMatrixD(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESTrp - ESIso')
title('Trp Mean Dark')
savefig('ShiMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg
%%
figure
imdata = log2(rdivide(VBEMatrixL(:,:,2),VBEMatrixL(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Light')
savefig('ShiVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = log2(rdivide(VBEMatrixD(:,:,2),VBEMatrixD(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Dark')
savefig('ShiVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = log2(rdivide(VBEMatrixL(:,:,3),VBEMatrixL(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Trp Variance Light')
savefig('ShiVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = log2(rdivide(VBEMatrixD(:,:,3),VBEMatrixD(:,:,1)));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Trp Variance Dark')
savefig('ShiVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

%%

medMatrix = abs(medianMatrix);
VBMatrix = abs(VBEMatrix);

figure
imdata = diff(medMatrix(:,:,[2 1]),[],3);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.3 0.3])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Shifts')
savefig('ShiMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = diff(medMatrix(:,:,[3 1]),[],3);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im2 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.3 0.3])
title(hColor,'ESTrp - ESIso')
title('Trp Mean Shifts')
savefig('TrpMedianMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = diff(VBMatrix(:,:,[2 1]),[],3);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Shifts')
savefig('ShiVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = diff(VBMatrix(:,:,[3 1]),[],3);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im4 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESTrp - ESIso')
title('Trp Variance Shifts')
savefig('TrpVBEMap.fig')
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(VBEMatrix(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(VBEMatrix(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
'CapSize',10)
yLab = ylabel('Variability Beyond\newline .   Expectation');
set(gca,'YTick',-0.4:0.2:0.6,'YLim',[-0.4 0.6],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
shg 
%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(medianMatrix(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(medianMatrix(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
'CapSize',10)
yLab = ylabel('Effect Size of Median');
set(gca,'YTick',-0.4:0.2:0.5,'YLim',[-0.5 0.5],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
shg

%% Make correlation heatmaps between all 5 metrics

idx = strcmp(array(:,4),'ISO');
blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);
mokidugway=interp1([1 128 129 256],[0 0.3 1; .9921 1 1; 1 .9921 .9921; 1 0 0],1:256);

YNames = {'Activity','Absolute Turn Bias', 'Turn Direction Switchiness','Wall Proximity',...
    'Turn Timing Clumpiness'};

data1 = data(idx,:);
error1 = error(idx,:);

genos = array(idx,3);
[a,b,c] = unique(genos);

means = nan(length(a),20);
for ii = 1:length(a)
means(ii,:) = nanmean(data1(c==ii,:));

data1(c==ii,:) = data1(c==ii,:)./means(ii,:);
end

data1(:,[12 17]) = abs(data1(:,[12 17])-1);

figure
corrs1 = corr(data1(:,16:20),'type','Spearman','rows','pairwise');
imagesc(corrs1)
colormap(mokidugway)
hCb = colorbar;
caxis([-1 1])
set(gca,'YTickLabel',YNames,'XTickLabel',YNames,'YTick',1:5,'XTickLabelRotation',25,'XTick',1:5,...
    'FontWeight','bold','FontSize',20)
title('Dark')
pbaspect([1 1 1])
ylabel(hCb,'\rho_{dark}','Rotation',270,'VerticalAlignment','cap')

figure
corrs2 = corr(data1(:,11:15),'rows','pairwise','type','Spearman');
imagesc(corrs2)
colormap(mokidugway)
hCb = colorbar;
caxis([-1 1])
set(gca,'YTickLabel',YNames,'XTickLabel',YNames,'YTick',1:5,'XTickLabelRotation',25,'XTick',1:5,...
    'FontWeight','bold','FontSize',20)
title('Light')
pbaspect([1 1 1])
ylabel(hCb,'\rho_{light}','Rotation',270,'VerticalAlignment','cap')

figure
imagesc(corrs2-corrs1)
colormap(mokidugway)
hCb = colorbar;
caxis([-0.5 0.5])
set(gca,'YTickLabel',YNames,'XTickLabel',YNames,'YTick',1:5,'XTickLabelRotation',25,'XTick',1:5,...
    'FontWeight','bold','FontSize',20)
title('Light - Dark')
pbaspect([1 1 1])
ylabel(hCb,'\rho_{light} - \rho_{dark}','Rotation',270,'VerticalAlignment','cap')
%% Make KDE for Delta Stats
ibecell = cell(0);
effectors = {'ISO','SHI','TRP'};
dates = cat(1,array{:,9});

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);


%  dataMatrix = nan(23,5,3);
%  errorMatrix = nan(23,5,3);
 
x = ones(size(array,1),1);%strcmp(dates(:,2),'02') & cell2mat(array(:,7)) == 4;
for ii = 1:size(importedData.data,1)
    %for kk = 1:4
    ngn = genonames(logical(importedData.data(ii,:)));
    neuron = ismember(upper(array(:,3)),ngn);
    
    idx = strcmp(array(:,4),'ISO') & neuron; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_ISO')%,'_tray',num2str(kk))
        ibescell{ii,1} = generateDiffPlot2(data(idx,:),error(idx,:));
        savefigs(name)
        close all
        dataMatrix(ii,:,1) = ibescell{ii,1}(:,1);
        errorMatrix(ii,:,1) = ibescell{ii,1}(:,2);
    end
    idx = strcmp(array(:,4),'SHI') & neuron; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_SHI')%,'_tray',num2str(kk))
        ibescell{ii,2} = generateDiffPlot2(data(idx,:),error(idx,:));
        savefigs(name)
        close all
        dataMatrix(ii,:,2) = ibescell{ii,2}(:,1);
        errorMatrix(ii,:,2) = ibescell{ii,2}(:,2);
    end
    idx = strcmp(array(:,4),'TRP') & neuron;%ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_TRP')%,'_tray',num2str(kk))
        ibescell{ii,3} = generateDiffPlot2(data(idx,:),error(idx,:));
        savefigs(name)
        close all
        dataMatrix(ii,:,3) = ibescell{ii,3}(:,1);
        errorMatrix(ii,:,3) = ibescell{ii,3}(:,2);
    end
  %  end
    
end

%dataMatrix(any(isnan(dataMatrix(:,:,2)) | isnan(dataMatrix(:,:,1)),2),:,:) = [];

%%

dataMatrix(dataMatrix==0) = nan;

figure
imdata = log2(dataMatrix(:,:,2)./dataMatrix(:,:,1));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
imagesc(imdata,'AlphaData',imAlpha)
set(gca,'YTickLabel',names,'YTick',1:size(dataMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0])
hColor = colorbar;
set(hColor,'Location','eastoutside','Color','black')
caxis([-0.8 0.8])
title(hColor,'Log2[shi/iso]')
title('Shibire Rank Order Shifts')
savefig('ShiRankMap.fig')
shg
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PRGn'))

figure
imdata = log2(dataMatrix(:,:,3)./dataMatrix(:,:,1));
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
imagesc(imdata,'AlphaData',imAlpha)
set(gca,'YTickLabel',names,'YTick',1:size(dataMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.8 0.8])
title(hColor,'Log2[trp/iso]')
title('Trp Rank Order Shifts')
savefig('TrpRankMap.fig')
shg
pbaspect([1 3.5 1])
colormap(gcf,brewermap(256,'PRGn'))

%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(dataMatrix(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(dataMatrix(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
'CapSize',10)
yLab = ylabel('LDM')
set(gca,'YTick',0:0.1:0.25,'YLim',[0 0.25],'FontSize',25,'XTick',[])

%% Make low > high shift figures
n = 6;
behave = 2;
[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(array(:,3)),ngn);
effectors = array(:,4);
generateShiftPlot(effectors(idx),data(idx,:),error(idx,:),behave)
title(names(n))
            xlabel('Proportion Pre-temp "Light-like"')


%% Generate a switchiness plot to estimate extinction coefficient

generateSwitchPlot(array(idx,:),lightDat);

          

%% Generate a transition triggered average
idx = strcmp(array(:,4),'ISO');
%figure
colors = {[0 0 0],[0.8 0.2 0],[0 0.6 0.2]}; 
alpha = 0.5;

figure
hold on

n = 15;
behave = 4;
light = 1;

f = @(p,x) p(1) + ((p(2)-p(1)) ./ ((1 + p(3)*exp(-p(4)*x)).^(1/p(5))));
h = fittype(@(p1,p2,p3,p4,p5,x) p1 + ((p2-p1) ./ ((1 + p3*exp(-p4*x)).^(1/p5))), 'options',s);

%f = @(p,x) p(1) + (p(2) - p(1)) ./ (1 + (x./p(3)).^p(4))%...
   % + p(5) + (p(6) - p(5)) ./ (1 + (x./p(7)).^p(8))

% f2 = @(p,x) p(1) + ((p(2)-p(1)) ./ ((p(3) + p(4)*exp(-p(5)*x)).^(1/p(6))))+p(7)*exp(-x);
%f = @(p,x) p(1)+((p(2)-p(1))./(1+(x./p(3)).^p(4)))+p(5)+((p(6)-p(5))./(1+(x./p(7)).^p(8)));
%f2 = @(p,x) p(1)*x + p(2);

% y = piecewise(x<0,f1,x>0,f2)

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));

idx = ismember(upper(array(:,3)),ngn) & strcmp(array(:,4),'ISO');
subArray = array(idx,:);
hError1 = generateTTAPlot(subArray,lightDat,behave,colors{1},light);
ind = hError1.XData > 0;
xVect = hError1.XData(1):0.1:hError1.XData(end);
a1 = fit(hError1.XData',hError1.YData',h);
hPlot1 = plot(xVect,a1(xVect),'Color',colors{1},'LineWidth',2);
hLine1 = line([a1(3), a1(3), -5, a1(3)],[0, f(a1,a1(3)), f(a1,a1(3)), f(a1,a1(3))],...
    'Color',[colors{1} alpha],'LineWidth',2,'LineStyle','-','LineJoin','round','LineWidth',0.5);


drawnow

idx = ismember(upper(array(:,3)),ngn) & strcmp(array(:,4),'SHI');
subArray = array(idx,:);
hError2 = generateTTAPlot(subArray,lightDat,behave,colors{2},light);
a2 = fit(hError2.XData',hError2.YData',h);
xVect = hError2.XData(1):0.1:hError2.XData(end);
hPlot2 = plot(xVect,a2(xVect),'Color',colors{2},'LineWidth',2);
hLine2 = line([a2(3), a2(3), -5, a2(3)],[0, f(a2,a2(3)), f(a2,a2(3)), f(a2,a2(3))],...
    'Color',[colors{2} alpha],'LineWidth',2,'LineStyle','-','LineJoin','round','LineWidth',0.5);


drawnow
% 
% idx = ismember(upper(array(:,3)),ngn) & strcmp(array(:,4),'TRP');
% if sum(idx)>0
% subArray = array(idx,:);
% hError3 = generateTTAPlot(subArray,lightDat,behave,colors{3},light);
% [a3,b] = sigm_fit(hError3.XData,hError3.YData,[],[0.08,0.08,0,0.5],0);
% xVect = linspace(hError3.XData(1),hError3.XData(end),100);
% hPlot3 = plot(xVect,f(a3,xVect),'Color',colors{3},'LineWidth',2);
% hLine3 = line([a3(3), a3(3), -20, a3(3)],[0, f(a3,a3(3)), f(a3,a3(3)), f(a3,a3(3))],...
%     'Color',[colors{2} alpha],'LineWidth',2,'LineStyle','-','LineJoin','round');
% drawnow
% end

title(names(n))
xlabel('Time relative to light change')
ylabel('Light-dependent modulation from light bias')
set(gca,'XLim',[-5 15],'TickDir','out','FontWeight','bold','FontSize',16,'YLim',[0 0.18])
pbaspect([1.5 1 1])
shg

%% Generate mutual information between light and turn direction as a function of time back
idx = strcmp(array(:,4),'ISO');
light = 0;
figure
%generateTurnTMIPlot(array(idx,:),lightDat,behave)
n = 5;
clf
hold on

ibecell = cell(0);
effectors = {'ISO','SHI','TRP'};
dates = cat(1,array{:,9});

[ug,~,ugi] = uniqueRowsCA(upper(array(:,3)));


%  dataMatrix = nan(23,5,3);
%  errorMatrix = nan(23,5,3);
 
x = ones(size(array,1),1);%strcmp(dates(:,2),'02') & cell2mat(array(:,7)) == 4;

    %for kk = 1:4
    ngn = genonames(logical(importedData.data(n,:)));
    neuron = ismember(upper(array(:,3)),ngn);
    idx = strcmp(array(:,4),'ISO') & neuron; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(n), '_ISO')%,'_tray',num2str(kk))
        hTurn1 = generateTurnTMIPlot(array(idx,:),lightDat,light,colors{1});
        val = (max(hTurn1.mainLine.YData)-min(hTurn1.mainLine.YData))/2+min(hTurn1.mainLine.YData);
        [~,ix] = min(abs(hTurn1.mainLine.YData(hTurn1.mainLine.XData <= 0) - val));
        line([hTurn1.mainLine.XData(1), hTurn1.mainLine.XData(ix), hTurn1.mainLine.XData(ix),...
            hTurn1.mainLine.XData(ix)],...
            [hTurn1.mainLine.YData(ix), hTurn1.mainLine.YData(ix), hTurn1.mainLine.YData(ix), 0],...
            'Color',colors{1})
        drawnow
    end
    idx = strcmp(array(:,4),'SHI') & neuron; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(n), '_SHI')%,'_tray',num2str(kk))
        hTurn2 = generateTurnTMIPlot(array(idx,:),lightDat,light,colors{2});
        val = (max(hTurn2.mainLine.YData)-min(hTurn2.mainLine.YData))/2+min(hTurn2.mainLine.YData);
        [~,ix] = min(abs(hTurn2.mainLine.YData(hTurn2.mainLine.XData <= 0) - val));
        line([hTurn2.mainLine.XData(1), hTurn2.mainLine.XData(ix), hTurn2.mainLine.XData(ix),...
            hTurn2.mainLine.XData(ix)],...
            [hTurn2.mainLine.YData(ix), hTurn2.mainLine.YData(ix), hTurn2.mainLine.YData(ix), 0],...
            'Color',colors{2})
        drawnow
    end
    idx = strcmp(array(:,4),'TRP') & neuron;%ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(n), '_TRP')%,'_tray',num2str(kk))
        hTurn3 = generateTurnTMIPlot(array(idx,:),lightDat,light,colors{3});
        val = (max(hTurn3.mainLine.YData)-min(hTurn3.mainLine.YData))/2+min(hTurn3.mainLine.YData);
        [~,ix] = min(abs(hTurn3.mainLine.YData(hTurn3.mainLine.XData <= 0) - val));
        line([hTurn3.mainLine.XData(1), hTurn3.mainLine.XData(ix), hTurn3.mainLine.XData(ix),...
            hTurn3.mainLine.XData(ix)],...
            [hTurn3.mainLine.YData(ix), hTurn3.mainLine.YData(ix), hTurn3.mainLine.YData(ix), 0],...
            'Color',colors{3})
        drawnow
    end
  %  end
  
  title(gca,names(n),'FontWeight','bold','FontSize',16)
  set(gca,'FontSize',12)
  xlabel(gca,'Time relative to turn')
  ylabel(gca,'Information (bits)')
    shg

%% Make LDM Bar plots

% import rubin annotation
filepath = '/Users/kyobikakaria/Desktop/Data/split_Rubin_Annotation.csv';
importedData = importdata(filepath);

names = importedData.textdata(2:end,1);
genonames = importedData.textdata(1,2:end);

shorterNames = {'P-FN12','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','PF-FGall','PF-FRub',...
    'PF-LCre','LPs-P',...
    'PBCap','PBInt','Ps-P1','Ps-P2'};

grps = [1 1 1 1 2 2 2 2 3 3 3 4 5 5 4 4 6 6 7 7 7 8];

names(1:length(shorterNames)) = shorterNames;

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);
[ue,~,uei] = unique(array(:,4));

screenCell = cell(length(ug),3);
neuronCell = cell(length(names),3);

VBEMean = cell(1,1);
VBEError = cell(1,1);
for jj = 1:size(ug,1)
    jj
    for kk = 1:size(ue,1)
        idx = ugi == jj & uei == kk;
        if sum(idx) > 0
        for ii = 1:5
            VBE = bootstrp(1000,@(x,y) IBEV2(x,y),data1(idx,[ii+10 ii+15]),...
                error1(idx,[ii+10 ii+15]));
            VBEMeanHigh{jj,kk}(ii) = nanmean(VBE);
            VBEErrorHigh{jj,kk}(ii) = nanstd(VBE);
            numFliesHigh(jj,kk) = sum(idx);
        end
        
        for ii = 1:5
            VBE = bootstrp(1000,@(x,y) IBEV2(x,y),data1(idx,[ii ii+5]),...
                error1(idx,[ii+10 ii+15]));
            VBEMeanLow{jj,kk}(ii) = nanmean(VBE);
            VBEErrorLow{jj,kk}(ii) = nanstd(VBE);
            numFliesLow(jj,kk) = sum(idx);
        end
        end
    end
end



%%
% genoPlots
clf

for ii = 1:size(VBEMeanHigh,1)
    ii
    subplot(8,7,ii)
    hold on
    dat = cat(1,VBEMeanHigh{ii,:});
    err = cat(1,VBEErrorHigh{ii,:});
    
    hBar = bar(dat');
    drawnow
    xData = repmat(1:size(dat,2),size(dat,1),1)+cat(1,hBar.XOffset);
    errorbar(xData',dat',err','LineStyle','none','CapSize',0,'LineWidth',2,'Color','black')
    

    title(ug(ii),'FontSize',10)
    set(gca,'XTickLabel',[])
    ylim([0 1])
end
shg



%% calculate deltas

for ii = 1:size(VBEMeanHigh,1)
    
    for jj = 1:2
        if ~isempty(VBEMeanHigh{ii,jj+1}) && ~isempty(VBEMeanHigh{ii,1})
            VBEMeanHighD{ii,jj} = (cell2mat(VBEMeanHigh(ii,jj+1)) - cell2mat(VBEMeanHigh(ii,1)))./...
                 cell2mat(VBEMeanHigh(ii,1));
             VBEErrorHighD{ii,jj} = sqrt((cell2mat(VBEErrorHigh(ii,jj+1)).^2 + ...
                 cell2mat(VBEErrorHigh(ii,1)).^2)./...
                 cell2mat(VBEMeanHigh(ii,1)));
        end
    end
    
end

%%
% neuronPlots

clf
neuronMeans = nan(16,5,2);
neuronErrors = nan(16,5,2);
for ii = 1:size(neuronCell,1)
    ngn = genonames(logical(importedData.data(ii,:)));
    [a,b,c] = intersect(ngn,upper(ug));
    if ~isempty(c)
    tempMeans = VBEMeanHighD(c,:);
    tempErrors = VBEErrorHighD(c,:);
    for jj = 1:size(tempMeans,2)
        if tempNums(:,jj) > 0
        neuronMeans(ii,:,jj) = nanmean(cat(1,tempMeans{:,jj}),1);
        e = cell2mat(tempErrors(:,jj));
        neuronErrors(ii,:,jj) = sqrt(sum(e.^2,1)./size(e,1));
        end
    end
    end
    subplot(8,2,ii)
    dat = squeeze(neuronMeans(ii,:,:));
    err = squeeze(neuronErrors(ii,:,:));
    
    hBar = bar(dat);
    drawnow
    hold on
    xData = repmat(1:size(dat,1),size(dat,2),1)+cat(1,hBar.XOffset);
    errorbar(xData,dat',err','LineStyle','none','CapSize',0,'LineWidth',0.5,'Color','black')
    
    set(gca,'XTickLabel',[],'YLim',[-0.6 0.3])
    title(names(ii))
    
    
    for kk = 1:size(dat,1)
    for ll = 1:size(dat,2)
        [~,p(kk,ll)] = ztest(dat(kk,ll),0,err(kk,ll))
        if p(kk,ll) < 0.001
            sym = '***';
        elseif p(kk,ll) < 0.01
            sym = '**';
        elseif p(kk,ll) < 0.05
            sym = '*';
        else
            sym = '';
        end
        
    end
end
    
    
end
%%

%%
    
ldm = bootstrp(1000,@(x) nanmean(abs(x(:,ii+10) - x(:,ii+15)))-nanmean(x(:,ii+10) - x(:,ii+15)),data1);
null = @(x,y) normrnd((x(:,ii+10)+x(:,ii+15))/2,sqrt(y(:,ii+10).^2+y(:,ii+15).^2)/2);
nulls = [null(data1,error1) null(data1,error1)];
ldmNull = bootstrp(1000,@(x,y) nanmean(abs(x(:,1) - x(:,2))) - nanmean(x(:,1) - x(:,2)),nulls);
ldmMean(:,ii) = nanmean(ldm)-nanmean(ldmNull);
ldmErr(:,ii) = nanstd(ldm);



%% 
for ii = 1:5
ldm = bootstrp(1000,@(x) nanmean(abs(x(:,ii+10) - x(:,ii+15))./x(:,ii+15)),data1);
null = @(x,y) normrnd((x(:,ii+10)+x(:,ii+15))/2,sqrt(y(:,ii+10).^2+y(:,ii+15).^2)/2);
ldmNull = bootstrp(1000,@(x,y) nanmean(abs(null(x,y) - null(x,y))./null(x,y)),data1,error1);
ldmMean(:,ii) = nanmean(ldm)-nanmean(ldmNull);
ldmErr(:,ii) = nanstd(ldm);
end

%%

ldmAll = bootstrp(1000, @(x)  nanmean(abs(mahal(x(:,11:15),x(:,11:15)) - ...
    mahal(x(:,16:20),x(:,16:20)))),data1)
nullAll = @(x,y) normrnd((x(:,11:15)+x(:,16:20))/2,sqrt(y(:,11:15).^2+y(:,16:20).^2)/2);
ldmNullAll = bootstrp(1000000,@(x,y) nanmean(abs(mahal(x,x) - ...
    mahal(y,y))),nullAll(data1,error1),nullAll(data1,error1));

idx = ~isnan(x) & ~isnan(y);

x = x(idx);
y = y(idx);

vals = 0.001:0.001:1
pd1 = pdf('Normal',vals,x,0.01);
pd2 = pdf('Normal',vals,x-y,0.01);

xy = floor([sum(pd1)' sum(pd2)'])


for ii = 1:size(pd1,1)
z(:,:,ii) = hist3([pd1(ii,:); pd2(ii,:)]');
end
%%
[a,~,c] = uniqueRowsCA(array(:,[3 4]));
data1 = cell2mat(array(:,15));
data2 = data1./nanmean(data1);
error1 = cell2mat(array(:,16));
error2 = error1./nanmean(data1);

for ii = 1:5
subplot(1,5,ii)
%scatter(data1(:,ii),data1(:,ii+5))
hist([data1(:,ii),data1(:,ii+5)])
end


f = @(x) corr(x(:,interleave2(11:15, 16:20)),'rows','pairwise');
bootstrp(100,@(x) f(normrnd(x{1},x{2})),{data2, error2})

corrmat = corr(data2(:,interleave2(11:15, 16:20)),'rows','pairwise');

blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);

imagesc(corrmat)

colormap(blanding)
caxis([-1 1])


shg
%% Generate simple figures

% this makes a coefficient of variation plot across all behavioral statistics
CV = generateCVPlot(array);

% this generates a correlation plot across all statistics
LDM = generateLDMPlot(array);


