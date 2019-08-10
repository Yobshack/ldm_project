%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 02.14.2018 - updated


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
lightDat = dlmread('lightSequenceScreen.txt');
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


%% Make LDM Bar plots

% import rubin annotation
filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/split_Gal4/split_Rubin_Annotation.csv';
importedData = importdata(filepath);

names = importedData.textdata(2:end,1);
genonames = importedData.textdata(1,2:end);


shorterNames = {'P-FN12','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','P-F-G','P-F-R',...
    'PF-LC','LPs-P',...
    'PBCap','Delta7','Ps-P1','Ps-P2'};

grps = [1 1 1 1 2 2 2 2 3 3 3 4 5 5 4 4 6 6 7 7 7 8];

names(1:length(shorterNames)) = shorterNames;

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);
[ue,~,uei] = unique(array(:,4));

screenCell = cell(length(ug),3);
neuronCell = cell(length(names),3);
%% Make large cell array to handle all data that I can query for groups

files = dir('*_processed.mat');
[allDataCell, cellColNames] = parseProcessedFiles(files);

%% Calculate behavioral metrics
% Takes a while to run this on a large data set
% order is light - low, dark - low, light - high, dark - high
% if ~isempty(dir('check1.mat'))
%     load('check1.mat')
% else
times = {[0 90]*60,[120 180]*60};
turnFilter = 100;
bins = 1;
allData = cell(0); allError = allData;
parfor hh = 1:size(allDataCell,1)
    hh
    [allData{hh}, allError{hh}] = generateBehavioralMetrics3(allDataCell(hh,:),lightDat,times,...
        bins,turnFilter);
    hh
end
%cellColNames(13:14) = {'LightLow','DarkLow','LightHigh','DarkHigh'};

%%
array = allDataCell;

array(:,13) = allData;
array(:,14) = allError;

% correct labeling for different genotypes
array(strcmp(array(:,3),'GMR-Gal4'),3) = cellstr('GMR');
array(strcmp(array(:,3),'Cry-39'),3) = cellstr('Cryptochrome');
array(strcmp(array(:,3),'cry-Gal4'),3) = cellstr('Cryptochrome');
array(strcmp(array(:,3),'PanR7-Gal4'),3) = cellstr('PanR7');
array(strcmp(array(:,3),'PanR8-Gal4'),3) = cellstr('PanR8');
array(strcmp(array(:,3),'Rh1-Gal4'),3) = cellstr('Rh1');
array(strcmp(array(:,3),'Rh3-Gal4'),3) = cellstr('Rh3');
array(strcmp(array(:,3),'Rh4-Gal4'),3) = cellstr('Rh4');
array(strcmp(array(:,3),'Rh5-Gal4'),3) = cellstr('Rh5');

array;
array(:,3) = upper(array(:,3));

data = cat(1,array{:,13});
error = cat(1,array{:,14});

idx = strcmp(array(:,4),'ISO');



%% Mean Shift plots
effectors = {'ISO','SHI','TRP'};
dates = cat(1,array{:,9});

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);


medianMatrix = nan(23,5,3,2);
VBEMatrix = nan(23,5,3,2);



subIdx = contains(genonames,'SS');
c = ismember(array(:,3),genonames(subIdx));
subArray = array(c,:);
subData = data(c,:);
subError = error(c,:);

medianCell = cell(0);
VBECell = cell(0);

%pca(subData(:,1:5))

x = ones(size(array,1),1);%strcmp(dates(:,2),'02') & cell2mat(array(:,7)) == 4;

isoKernels = nan(length(names),1001,2);
shiKernels = isoKernels;
trpKernels = isoKernels;

for ii = 1:length(names)
    %     %for kk = 1:4
        ngn = genonames(logical(importedData.data(ii,:)));
        
%     neuron = ugi == ii;
    %idx = 1:size(data,1);
    temp = nan(4,5,length(ngn),3);
    temp2 = nan(8,5,length(ngn),1);
    for mm = 1:length(ngn)
    neuron = ismember(upper(subArray(:,3)),ngn(mm));
    idx = strcmp(subArray(:,4),'ISO') & neuron; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_Distributions_ISO')%,'_tray',num2str(kk))
        
        [vars,vars2] = generateDistPlots2(subData(idx,:),subError(idx,:));
        temp(:,:,mm,1) = vars;
        temp2(:,:,mm,1) = vars2;
    end
    idx = strcmp(subArray(:,4),'SHI') & neuron; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_Distributions_SHI')%,'_tray',num2str(kk))
        [vars,vars2] =  generateDistPlots2(subData(idx,:),subError(idx,:));
        temp(:,:,mm,2) = vars;

    end
    idx = strcmp(subArray(:,4),'TRP') & neuron;%ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_Distributions_TRP')%,'_tray',num2str(kk))
        [vars,vars2] = generateDistPlots2(subData(idx,:),subError(idx,:));
        
        temp(:,:,mm,3) = vars;
        
    end
    
    
    
    end
    
    f = @(x,y) (x-y)./y;
    
    x1 = nanmean(temp2(3:4,:,:),3);
    medianCell{ii,1} = x1;
    medianMatrix(ii,:,1,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    x1 = nanmean(temp2(7:8,:,:),3);
    VBECell{ii,1} = x1;
    VBEMatrix(ii,:,1,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    
    x1 = f(temp(1:2,:,:,2),temp(1:2,:,:,1));
    medianCell{ii,2} = x1;
    medianMatrix(ii,:,2,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    x1 = f(temp(3:4,:,:,2),temp(3:4,:,:,1));
    VBECell{ii,2} = x1;
    VBEMatrix(ii,:,2,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    
    x1 = f(temp(1:2,:,:,3),temp(1:2,:,:,1));
    medianCell{ii,3} = x1;
    medianMatrix(ii,:,3,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    x1 = f(temp(3:4,:,:,3),temp(3:4,:,:,1));
    VBECell{ii,3} = x1;
    VBEMatrix(ii,:,3,:) = reshape(nanmean(x1,3)',[1,5,1,2]);
    %  end
end

%% Make violin plots for a set of lines

figure
hold on
bins = 0:0.001:1;
line = 20
ug(line)
for jj = line
    for kk = 1:4
        for ii = 1:2
            if ii == 1
                color = c1;
            else
                color = c2;
            end
            if kk == 1
                y1 = isoKernels(jj,:,ii);
            elseif kk == 2
                y1 = shiKernels(jj,:,ii);
            elseif kk == 3
                y1 = trpKernels(jj,:,ii);
            elseif kk == 4
                y1 = shiKernels(jj,:,ii) - isoKernels(jj,:,ii);
            end
            y1 = y1./max(y1)*0.3;
            
            idx = y1 > 0.01;
            
            offset = kk*2-2+ii;
            hArea1 = fill([offset-y1(idx) fliplr(offset+y1(idx))],...
                [bins(idx) fliplr(bins(idx))],color);
        end
    end
end
%hArea2 = area(y1,bins,'LineWidth',2,'FaceAlpha',0.5,'FaceColor',c2);

set(gca,'YLim',[0.35 0.65])

%%
%         ngn = genonames(logical(importedData.data(ii,:)));
%         neuron = ismember(upper(array(:,3)),ngn);

figure
imdata = [log2(rdivide(medianMatrixL(:,:,2),medianMatrixL(:,:,1))),...
    log2(rdivide(medianMatrixD(:,:,2),medianMatrixD(:,:,1)))];
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',ug,'YTick',1:size(medianMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Light')
savefig('ShiMedianMap.fig')
pbaspect([1 16 1])
colormap(gcf,brewermap(256,'RdBu'))
shg
%%

figure
imdata = medianMatrix(:,:,1,:);
hBar = barwitherr(squeeze(nanstd(imdata,1)),squeeze(nanmean(imdata,1)));
set(gca,'XTick')

figure
imdata = VBEMatrix(:,:,1,:);
hBar = barwitherr(squeeze(nanstd(imdata,1)),squeeze(nanmean(imdata,1)));
set(gca,'XTick')

%%
figure
imdata = medianMatrix(:,:,2,1);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medianMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Light')
%savefig('ShiMedianMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = medianMatrix(:,:,2,2);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medianMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESTrp - ESIso')
title('Shibire Mean Dark')
%savefig('ShiMedianMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata =  medianMatrix(:,:,3,1);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medianMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESTrp - ESIso')
title('Trp Mean Light')
savefig('ShiMedianMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'RdBu'))
shg

figure
imdata = medianMatrix(:,:,3,2);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im1 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'blues'))
set(gca,'YTickLabel',names,'YTick',1:size(medianMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-0.5 0.5])
title(hColor,'ESTrp - ESIso')
title('Trp Mean Dark')
savefig('ShiMedianMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'RdBu'))
shg
%%
figure
imdata = [log2(rdivide(VBEMatrixL(:,:,2),VBEMatrixL(:,:,1))),...
    log2(rdivide(VBEMatrixD(:,:,2),VBEMatrixD(:,:,1)))];
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',ug,'YTick',1:length(ug),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Light')
savefig('ShiVBEMap.fig')
pbaspect([1 16 1])
colormap(gcf,brewermap(256,'PuOr'))
shg
%%

figure
imdata = VBEMatrix(:,:,2,1);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBEMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Light')
savefig('ShiVBEMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'PuOr'))
shg


figure
imdata = VBEMatrix(:,:,2,2);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBEMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Shibire Variance Dark')
savefig('ShiVBEMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = VBEMatrix(:,:,3,1);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBEMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Trp Variance Light')
savefig('ShiVBEMap.fig')
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'PuOr'))
shg

figure
imdata = VBEMatrix(:,:,3,2);
imAlpha = ones(size(imdata));
imAlpha(isnan(imdata))=0;
im3 = imagesc(imdata,'AlphaData',imAlpha);
colormap(brewermap(12,'reds'))
set(gca,'YTickLabel',names,'YTick',1:size(VBEMatrix,1),'MinorGridLineStyle','none','Box','off',...
    'TickLength',[0 0],'Color','black')
hColor = colorbar;
set(hColor,'Location','eastoutside')
caxis([-1 1])
title(hColor,'ESShi - ESIso')
title('Trp Variance Dark')
savefig('ShiVBEMap.fig')
pbaspect([1 3 1])
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

x = nanmean(VBEMatrixL(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(VBEMatrixL(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
    'CapSize',10)
yLab = ylabel('Variability Beyond\newline .   Expectation');
set(gca,'YTick',-0.4:0.2:0.6,'YLim',[-0.4 0.6],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
title('VBEMatrixL')
shg

%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(VBEMatrixD(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(VBEMatrixD(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
    'CapSize',10)
yLab = ylabel('Variability Beyond\newline .   Expectation');
set(gca,'YTick',-0.4:0.2:0.6,'YLim',[-0.4 0.6],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
title('VBEMatrixD')
shg

%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(medianMatrixL(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(medianMatrixL(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
    'CapSize',10)
yLab = ylabel('Variability Beyond\newline .   Expectation');
set(gca,'YTick',-0.4:0.2:1,'YLim',[-0.4 1],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
title('medianMatrixL')
shg

%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(medianMatrixD(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(medianMatrixD(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
    'CapSize',10)
yLab = ylabel('Variability Beyond\newline .   Expectation');
set(gca,'YTick',-0.4:0.2:1,'YLim',[-0.4 1],'FontSize',25,'XTick',[],'XLim',[0 6])
pbaspect([1 1 1])
title('medianMatrixD')
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

data1(:,[2 7]) = abs(data1(:,[2 7])-1);

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
ibescell = cell(0);
effectors = {'ISO','SHI','TRP'};
dates = cat(1,array{:,9});

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);


 dataMatrix = nan(23,5,2);
 errorMatrix = nan(23,5,2);

% dataMatrix = nan(23,1,3);
% errorMatrix = nan(23,1,3);
% f1 = figure;
% f2 = figure;


subIdx = contains(genonames,'SS');
c = ismember(array(:,3),genonames(subIdx));
subArray = array(c,:);
subData = data(c,:);
subError = error(c,:);

x = ones(size(array,1),1);%strcmp(dates(:,2),'02') & cell2mat(array(:,7)) == 4;
for ii = 1:16
    %for kk = 1:4
    ngn = genonames(logical(importedData.data(ii,:)));
    
    tempMatrix1 = nan(length(ngn),5,3);
    tempMatrix2 = tempMatrix1;
    for mm = 1:length(ngn)
        neuron = ismember(upper(subArray(:,3)),ngn(mm));
   % neuron = ugi == ii;
    idx = strcmp(subArray(:,4),'ISO') & neuron; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_ISO')%,'_tray',num2str(kk))
        ibescell{ii,1} = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,1) = ibescell{ii,1}(:,1);
        tempMatrix2(mm,:,1) = ibescell{ii,1}(:,2);
    end
    idx = strcmp(subArray(:,4),'SHI') & neuron; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_SHI')%,'_tray',num2str(kk))
        ibescell{ii,2} = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,2) = ibescell{ii,2}(:,1);
        tempMatrix2(mm,:,2) = ibescell{ii,2}(:,2);
    end
    idx = strcmp(subArray(:,4),'TRP') & neuron;%ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_TRP')%,'_tray',num2str(kk))
        ibescell{ii,3} = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,3) = ibescell{ii,3}(:,1);
        tempMatrix2(mm,:,3) = ibescell{ii,3}(:,2);
    end
    %  end
    end
      f = @(x,y) (x-y)./y;
    out(ii,:) = nanmean(tempMatrix1(:,:,1));
    dataMatrix(ii,:,1) = nanmean(f(tempMatrix1(:,:,2),tempMatrix1(:,:,1)),1);
    dataMatrix(ii,:,2) = nanmean(f(tempMatrix1(:,:,3),tempMatrix1(:,:,1)),1);
    

end

%dataMatrix(any(isnan(dataMatrix(:,:,2)) | isnan(dataMatrix(:,:,1)),2),:,:) = [];



%%

dataMatrix(dataMatrix==0) = nan;

figure
imdata = dataMatrix(:,:,1);%(dataMatrix(:,:,2)-dataMatrix(:,:,1))./dataMatrix(:,:,1);
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
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'PRGn'))

figure
imdata = dataMatrix(:,:,2);%(dataMatrix(:,:,3)-dataMatrix(:,:,1))./dataMatrix(:,:,1);
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
pbaspect([1 3 1])
colormap(gcf,brewermap(256,'PRGn'))

%% Make KDE for Delta Stats

% cutoffs = [0:0.05:0.5];
% global cutoff
% for kk = 1:length(cutoffs)
ibescell = cell(0);
effectors = {'ISO','SHI','TRP'};
dates = cat(1,array{:,9});

[ug,~,ugi] = uniqueRowsCA([upper(array(:,3))]);


%  dataMatrix = nan(23,5,2);
%  errorMatrix = nan(23,5,2);

% dataMatrix = nan(23,1,3);
% errorMatrix = nan(23,1,3);
% f1 = figure;
% f2 = figure;


subIdx = contains(genonames,'SS');
%subIdx = ones(length(genonames),1);
c = ismember(array(:,3),genonames(subIdx));
subArray = array(c,:);
subData = data(c,:);
subError = error(c,:);


testNeuron = 7;
ii = testNeuron
x = ones(size(array,1),1);%strcmp(dates(:,2),'02') & cell2mat(array(:,7)) == 4;

%    cutoff = cutoffs(kk);
    %for kk = 1:4
  for ii = 1:16
    ngn = genonames(logical(importedData.data(ii,:)));
    %tray = subArray



    tempMatrix1 = nan(length(ngn),5,3);
    tempMatrix2 = tempMatrix1;
    nMatrix = tempMatrix1;
    for mm = 1:length(ngn)
        neuron = ismember(upper(subArray(:,3)),ngn(mm));
   % neuron = ugi == ii;
    idx = strcmp(subArray(:,4),'ISO') & neuron; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_ISO')%,'_tray',num2str(kk))
        [temp,temp2] = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,1) = temp(1);
        tempMatrix2(mm,:,1) = temp(2);
        nMatrix(mm,:,1) = temp2;
    end
    idx = strcmp(subArray(:,4),'SHI') & neuron; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_SHI')%,'_tray',num2str(kk))
        [temp,temp2] = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,2) = temp(1);
        tempMatrix2(mm,:,2) = temp(2);
        nMatrix(mm,:,2) = temp2;
    end
    idx = strcmp(subArray(:,4),'TRP') & neuron;%ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        name = strcat(names(ii), '_TRP')%,'_tray',num2str(kk))
        [temp,temp2] = generateDiffPlot3(subData(idx,:),subError(idx,:));
%         savefigs(name)
%         close all
        tempMatrix1(mm,:,3) = temp(1);
        tempMatrix2(mm,:,3) = temp(2);
        nMatrix(mm,:,3) = temp2;
    end
    %  end
        
    if tempMatrix2(mm,2,3) > 0.25
        tempMatrix1(mm,:,3) = nan;
        tempMatrix2(mm,:,3) = nan;
    end
    end
    
    f = @(x,y) (x-y)./y;
    out(ii,:) = nanmean(tempMatrix1(:,:,1),1);
    dataMatrix(ii,:,1) = nanmean(f(tempMatrix1(:,:,2),tempMatrix1(:,:,1)),1);
    f2 = @(x,y) sqrt(nanmean((x.^2+y.^2)./tempMatrix1(:,:,1).^2))...
        ./sqrt(sum(~isnan(tempMatrix1(:,2,2))));
    errorMatrix(ii,:,1) = f2(tempMatrix2(:,:,1),tempMatrix2(:,:,2));
    dataMatrix(ii,:,2) = nanmean(f(tempMatrix1(:,:,3),tempMatrix1(:,:,1)),1);
    errorMatrix(ii,:,2) = f2(tempMatrix2(:,:,1),tempMatrix2(:,:,2));
    
    ind = all(~isnan(tempMatrix1(:,2,:)),3);
    tempAll{ii} = tempMatrix1(ind,2,:);
    tempAllE{ii} = tempMatrix2(ind,2,:);
    nameAll{ii} = ngn(ind);

    
    clf
    c1 = [0 0.6 0]
    
%     hold on
%     hBar = bar(squeeze(nanmean(tempMatrix1(:,2,:),1)),'FaceColor',c1,'BarWidth',0.5,...
%         'LineWidth',3,'EdgeColor',[0 0 0],'EdgeAlpha',1)
       ind = ~isnan(tempMatrix1(:,2,1));
       f = @(x) (x-0.15):(0.3/(sum(ind)-1)):(x+0.15);
       x = squeeze(tempMatrix1(ind,2,1:2));
       y =  squeeze(tempMatrix2(ind,2,1:2));
       n = squeeze(nMatrix(ind,2,1:2));
       
       y(any(isnan(x),2),:) = nan;
       x(any(isnan(x),2),:) = nan;
      
     hError = errorbar([f(1);f(2)],x',...
       y','Marker','o','LineStyle','--','CapSize',0,...
        'MarkerSize',20,'MarkerFaceColor','auto','LineWidth',2);
    colormap('gray')
    effs = {'Gal4/+','Gal4/Shits'};
    set(gca,'XLim',[0.5 2.5],'YLim',[0.5 1],'XTick',[1:3],'XTickLabel',effs,'FontName','arial',...
        'FontWeight','bold','FontSize',28,'Color','none','YTick',[0:0.1:1])
    set(gcf,'Renderer','painters')
    ylabel('Normalized LDM')
    legend(ngn(ind))
    legend('boxoff')
    pbaspect([0.5 1 1])
     name = strcat(names{ii},'_ldms');
        title(name,'Interpreter','none')

    
    pvals = cell(0);
    if size(x,2) == 1
        x = x';
        y = y';
    end
    for jj = 1:size(x,1)
    pvals{jj} = num2str(flyVacPersMetricPValue([x(jj,1),y(jj,1),x(jj,2),y(jj,2)]));
    end
    
    text(1,0.6,strcat({pvals{:},num2str(n)}))
    drawnow
    shg
   
    savefig(gcf,name)
    print(name,'-dpdf','-fillpage')
    
    end

%   end
%% Make plot of data set to look at multiple neurons at once
clf
indxs = [1:8 12 15:16 14 9:11 19];
c = [0 0 0.4];
line([0 length(indxs)+1],[0 0],'Color','red','LineWidth',2);
hold on
errorbar(dataMatrix(indxs,2,1),errorMatrix(indxs,2,1),'Marker','o','MarkerSize',15,...
    'MarkerFaceColor','auto','CapSize',0,'LineWidth',3,'MarkerEdgeColor','auto','Color',c,...
    'LineStyle','none')
set(gca,'XLim',[0 length(indxs)+1],'XTickLabel',names(indxs),'XTickLabelRotation',45,'XTick',1:length(indxs),...
    'FontSize',16,'FontName','arial','YTickLabel',[-30:10:10],'FontWeight',...
    'bold','YTick',[-0.3:0.1:0.1])
ylabel('Percent LDM Change')
name = 'input_LDM';
savefig(gcf,name)
pbaspect([1.5 1 1])
print(name,'-dpdf','-fillpage')

shg

%% Analyze lines for certain neurons as a sliding window

neuronClass = 11;
importedData.textdata(neuronClass+1,1)

timeBin = 3600;
timeVect = 0:1800:10800;

turnfilter = 100;

effs = {'ISO','SHI','TRP'};

ngn = genonames(logical(importedData.data(neuronClass,:)));
ldm = nan(length(timeVect),length(effs),length(ngn));
ldmE = ldm;
nf = ldm;
covM = ldm;
products = ldm;


for mm = 1:length(ngn)
neuron = ismember(upper(array(:,3)),ngn(mm)); 

subArray = array(neuron,:);

tb = nan(size(subArray,1),length(timeVect),2);
tbe = tb;




for kk = 1:3
    idx = strcmp(subArray(:,4),effs(kk));
    if sum(idx) < 30
        continue
    end
    subArray2 = subArray(idx,:);
    for nn = 1:size(subArray2,1)
            times = subArray2{nn,2};
            turns = subArray2{nn,1};
            lights = subArray2{nn,11};

            
            [tRel tRelType] = calculateRelativeTime(times,lightDat);
            ind = tRel > 5 | tRel < -2;
            subArray2{nn,2} = times(ind);
            subArray2{nn,1} = turns(ind);
            subArray2{nn,11} = lights(ind);
    end
            
    for jj = 1:size(tb,2)
        for ii = 1:size(subArray2,1)
            times = subArray2{ii,2};
            turns = subArray2{ii,1};
            lights = subArray2{ii,11};

 
            tIdxL = (times > timeVect(jj)) & (times <= timeVect(jj) + timeBin) & lights;
            tempL = turnbias(turns(tIdxL));
            tIdxD = (times > timeVect(jj)) & (times <= timeVect(jj) + timeBin) & ~lights;
            tempD = turnbias(turns(tIdxD));
            
            if sum(tIdxL) > turnfilter && sum(tIdxD) > turnfilter 
            tb(ii,jj,:) = cat(3,tempL(1), tempD(1));
            tbe(ii,jj,:) = cat(3,tempL(2), tempD(2));

            end
        end
%         

        tempBoots = bootstrp(1000,@(x,y) IBEV5(x,y),cat(2,tb(:,jj,1),tb(:,jj,2)),...
            cat(2,tbe(:,jj,1),tbe(:,jj,2)));
%         tempBoots = bootstrp(1000,@(x,y) mad(x-y),tb(:,jj,1),tb(:,jj,2));
        ldm(jj,kk,mm) = nanmean(tempBoots);
        ldmE(jj,kk,mm) = nanstd(tempBoots);
        nf(jj,kk,mm) = sum(~isnan(tb(:,jj,1)));
        if nf(jj,kk,mm) > 0
            temp = cov(ldm(jj,kk,mm),ldm(1,kk,mm),'omitrows');
            covM(jj,kk,mm) = temp(1,2);
            products(jj,kk,mm) = nanmean(tb(:,jj,1).*tb(:,jj,2));
        end
        
        jj
        
    end
    

    
end
end
%%
clf
a = ldm;
temp = (a(:,2:3,:) - a(:,1,:))./a(:,1,:);
ldm2 = nanmean(temp,3);
b = abs(a).*sqrt( (ldmE./ldm).^2+(ldmE(1,:,:)./ldm(1,:,:)).^2 - ...
    2*covM./(ldm.*ldm(1,:,:)) );
b2 = (b(:,2:3,:).^2 + b(:,1,:).^2)./a(:,1,:).^2;
ldmE2 = sqrt(nanmean(b2,3));
c1 = [1 0 0]
c2 = [0 0 1]

timeBins = timeVect+0.5*timeBin;
for hh = 1:2
    
        switch hh
            case 1
                c = c1;
                varLeg = 'Shi';
            case 2
                c = c2;
                varLeg = 'Trp';
        end

    subplot(1,2,hh)
hold on
plot(repmat(timeBins/60,2,1)',ldm2(:,hh),'LineWidth',5,'Color',c)

for nn = 1:size(a,3)

plot(timeBins/60,temp(:,hh,nn),'Color',c)


text(timeBins(end)/60,temp(end,hh,nn),ngn(nn))

end

line([0 180],[0 0],'Color',[0 0 0]);
legend(varLeg)
title(names(neuronClass))
set(gca,'YLim',[-1 1])
end
shg
%% Make little bar plots to fit on top of heat maps
clf
hold on

x = nanmean(dataMatrix(:,:,1));
hBar1 = bar(x,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
errorbar(1:5,x,nanstd(dataMatrix(:,:,1)),'LineStyle','none','LineWidth',5,'Color','black',...
    'CapSize',10)
yLab = ylabel('LDM');
set(gca,'YTick',0:0.1:0.25,'YLim',[0 0.25],'FontSize',25,'XTick',[])

%%

subIdx = contains(genonames,'');
%subIdx = ones(length(genonames),1);
c = ismember(array(:,3),genonames(subIdx));
subArray = array(c,:);
subData = data(c,:);
subError = error(c,:);


%% Make low > high shift figures
n = 1;
%shiftMat = nan(23,4);
for n = 1:23
    clf
%subplot(4,6,n)
behave = 2;
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(subArray(:,3)),ngn);
% split = 'SS02296'
%idx = ismember(array(:,3),split) 
effectors = subArray(:,4);
name = strcat(names{n},'_shifts');
temp = generateShiftPlot(effectors(idx),subData(idx,:),subError(idx,:),behave);
shiftMat(n,:) = [diff(temp(:,1:2),[],2)' sqrt(sum(temp(:,3:4).^2,2))']
title(names(n))
ylabel({'PSI','Proportion Light-Like'})
drawnow
% savefig(gcf,name)
% print(name,'-dpdf','-fillpage')
end


%% Make plot of data set to look at multiple neurons at once
clf
indxs = [1:8 12 15:16 14 9:11 19];
c = [0 0 0.4];
mokidugway=interp1([1 128 129 256],[0 1 1; .9921 1 1; 1 .9921 .9921; 1 0 0],1:256);
colormap(mokidugway)
figure
y = shiftMat(indxs,1:2);
e = shiftMat(indxs,3:4);

% hBar = bar(y);
% hBar(1).FaceColor = [22/255 231/255 207/255];
% hBar(2).FaceColor = [0/255 118/255 186/255];
%     drawnow
%     hold on
%     for ii = 1:size(y,2)
%         hError = errorbar(hBar(ii).XData+hBar(ii).XOffset,y(:,ii),e(:,ii),...
%             'LineStyle','none','Color','black','CapSize',0,'LineWidth',1);
%     end
c = [0/255 118/255 186/255];;
errorbar(y(:,2),e(:,2),'Marker','o','MarkerSize',10,...
    'MarkerFaceColor','auto','CapSize',0,'LineWidth',3,'MarkerEdgeColor','auto','Color',c,...
    'LineStyle','none')
set(gca,'XLim',[0 length(indxs)+1],'XTickLabel',names(indxs),'XTickLabelRotation',45,'XTick',1:length(indxs),...
    'FontSize',16,'FontName','arial','FontWeight',...
    'bold','YLim',[-0.5 0.5],'YTick',-0.5:0.25:0.5)
hold on
line([0 length(indxs)+1],[0 0],'Color','red','LineWidth',2);
ylabel('PSI Shift')
name = 'PSIShift_Dark';
savefig(gcf,name)
pbaspect([1.5 1 1])
print(name,'-dpdf','-fillpage')


shg

%% Make low > high shift figures
n = 1;

behave = 2;
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(subArray(:,3)),ngn);
for ii = 1:length(ngn)
    clf
split = ug(ii)
split = 'SS02252'
idx = ismember(subArray(:,3),split);
effectors = subArray(:,4);
generateShiftPlot(effectors(idx),subData(idx,:),subError(idx,:),behave)
title(split)
ylabel('Pre-Temp Similarity Index')
xlabel('High Temp Condition')
%xlabel('Proportion Pre-temp "Light-like"')
savefig(gcf,split{:})
print(split{:},'-dpdf','-fillpage')
end

%%
n = 20;
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(subArray(:,3)),ngn) & strcmp(effectors,'SHI');
data = subData(idx,[2 7 12 17])
[a,b] = sort(data);
[b1, i1] = sort(b);
[~,temp] = min(abs(i1(:,1) - i1(:,3:4)),[],2);
lightSHI(:,1) = temp == 1;

[~,temp] = min(abs(i1(:,2) - i1(:,3:4)),[],2);
lightSHI(:,2) = temp == 1;
   
ind = find(diff(lightSHI,[],2)<0);

Names = {'B','C','D','E','F'};

subplot(1,4,1)
hB1 = barh(data(ind,1)-nanmean(dataLow(:,1)),'BarWidth',0.5)
subplot(1,4,2)
hB2 = barh(data(ind,2)-nanmean(data(:,2)),'BarWidth',0.5)
subplot(1,4,3)
hB3 = barh(data(ind,3)-nanmean(data(:,3)),'BarWidth',0.5)
subplot(1,4,4)
hB4 = barh(data(ind,4)-nanmean(data(:,4)),'BarWidth',0.5)
set([hB1.Parent,hB2.Parent,hB3.Parent,hB4.Parent],'XLim',[-0.35 0.35],...
    'PlotBoxAspectRatio',[0.5 1 1],'YTickLabel',Names)
ylabel([hB1.Parent],'Individual Flies')

shg

%% Make shift scatter plots
n = 20;
clf
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
ngn = 'SS02252'
idx = ismember(upper(subArray(:,3)),ngn) & ~strcmp(subArray(:,4),'TRP') &...
    all(~isnan(subData(:,2:5:end)),2) & all(subData(:,1:5:end) > 0.01,2);
data1 = subData(idx,2:5:end);


f = @(x) (x-nanmean(x,1))./nanstd(x)
%data1 = f(data1);
effs = subArray(idx,4);
[a,b,c] = unique(effs);
nyidalur=interp1([1 128 129 256],[1 0 1; .2559 .248 .2559; .248 .2559 .248; 0 1 0],1:256);
blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);

colormap([0 0 0;1 0 0])

cs = {[0 0 0],[1 0 0]}

f = @(m,x) m(1)*x+m(2); 
comps = {[3,1], [3,2], [4,2], [4,1]}
for ii = 1:4
    %subplot(2,2,ii)
    hold on
    x = data1(:,comps{ii}(1));
    y = data1(:,comps{ii}(2));
%     scatter(x,y,'CData',c,'MarkerFaceColor','flat',...
%         'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0,'SizeData',repmat(15,size(x,1),1))
%     [Z,A,B,ALPHA]  = fitellipse([x(c==1,:),y(c==1,:)],'linear');
%     hE1 = plotellipse(Z,A,B,ALPHA);
%     [Z,A,B,ALPHA]  = fitellipse([x(c==2,:),y(c==2,:)],'linear');
%     hE2 = plotellipse(Z,A,B,ALPHA);
% 
%     fit_ellipse(x(c==1,:),y(c==1,:),gca);
%     fit_ellipse(x(c==2,:),y(c==2,:),gca);

    
    cor(ii,1) = corr(x(c==1,:),y(c==1,:),'type','Spearman');
    cor(ii,2) = corr(x(c==2,:),y(c==2,:),'type','Spearman');
    
    corE(ii,1) = nanstd(bootstrp(1000,@(x,y) corr(x,y,'type','Spearman'),...
        x(c==1,:),y(c==1,:)));
    corE(ii,2) = nanstd(bootstrp(1000,@(x,y) corr(x,y,'type','Spearman'),...
        x(c==1,:),y(c==1,:)));
%     m1 = polyfit(x(c==1,:),y(c==1,:),1);
%     m2 = polyfit(x(c==2,:),y(c==2,:),1);
%     
% 
%     
%     xVals = [min(x),max(x)];
%     line(xVals,f(m1,xVals),'color','black');
%     line(xVals,f(m2,xVals),'color','red');
%     
%     set(gca,'XLim',[0 1],'YLim',[0 1],'FontName','arial')
%     pbaspect([1 1 1])
%     text(0.6, 0.1, [strcat({'r = '},num2str(round(cor1,2))),strcat({'r = '},num2str(round(cor2,2)))])
%     

end

y = cor;
e = corE;
hBar = bar(y);
drawnow

for ii = 1:size(y,2)
    hError = errorbar(hBar(ii).XData+hBar(ii).XOffset,y(:,ii),e(:,ii),...
        'LineStyle','none','Color','black','CapSize',0,'LineWidth',3);
end

name = strcat('PSI_Scatter');

set(gca,'XTickLabel',{'LL','LD','DD','DL'},'XTick',1:4,'FontName','arial')
print(name,'-dpdf','-bestfit')

shg

%% Generate Fig. 2A (Sliding Window)
n = [20];
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(any(logical(importedData.data(n,:)),1));
idx = ismember(upper(subArray(:,3)),ngn) & all(subData(:,[1 6]) > 0.03,2);
subArray2 = subArray(idx,:);

names(n)
timebin = 30;
dEnd = [0:5:(180-timebin)]*60;
timebin = timebin*60;
timeMat = nan(size(subArray2,1),length(dEnd));
x = nan(size(subArray2,1), 2);
y = x;
nullX = x;
nullY = x;
ldmIso = nan(length(dEnd)-1,2); ldmShi = ldmIso;

for jj = 1:length(dEnd)-1
    for ii = 1:size(subArray2,1)
        turns = subArray2{ii,1};
        times = subArray2{ii,2};
        lights = subArray2{ii,11};
        
        idx = times >= dEnd(jj) & times < dEnd(jj) + timebin;
        if sum(idx) > turnfilter
        x(ii,:) = turnbias(turns(idx & lights));
        y(ii,:) = turnbias(turns(idx & ~lights));
        
        if jj == 1
            nullX(ii,:) = x(ii,:);
            nullY(ii,:) = y(ii,:);
        end
        end
    end
    ind= strcmp(subArray2(:,4),'ISO');
    dat = [nullX(ind,1),nullY(ind,1),x(ind,1), y(ind,1)];
    err = [nullX(ind,2),nullY(ind,2),x(ind,2), y(ind,2)];
    temp1 = IBEV12(dat,err);
    temp2 = bootstrp(1000,@(x,y) IBEV12(x,y),dat,...
        err);
    ldmIso(jj,:) = [temp1 nanstd(temp2)];
    n = sum(all(~isnan(dat),2));
    
    ind= strcmp(subArray2(:,4),'SHI');
    dat = [nullX(ind,1),nullY(ind,1),x(ind,1), y(ind,1)];
    err = [nullX(ind,2),nullY(ind,2),x(ind,2), y(ind,2)];
    temp1 = IBEV12(dat,err);
    [temp2] = bootstrp(1000,@(x,y) IBEV12(x,y),dat,...
        err);
    ldmShi(jj,:) = [temp1 nanstd(temp2)];
    n = sum(all(~isnan(dat),2))
    
end

clf
xVect = (dEnd(1:(end-1))+timebin/2)/60;
timeI = xVect <= 60;
hPlot = plot(xVect,[ldmIso(:,1)./nanmean(ldmIso(timeI,1)),ldmShi(:,1)./...
    nanmean(ldmShi(timeI,1))],'LineWidth',3);

hPatch1 = patch([xVect fliplr(xVect)],[hPlot(1).YData'+ldmIso(:,2);...
    flipud(hPlot(1).YData'-ldmIso(:,2))],...
    hPlot(1).Color);
hPatch2 = patch([xVect fliplr(xVect)],[hPlot(2).YData'+ldmShi(:,2);...
    flipud(hPlot(2).YData'-ldmShi(:,2))],...
    hPlot(2).Color);

set([hPatch1,hPatch2],'FaceAlpha',0.25,'EdgeAlpha',0)
set(gca,'XLim',[0 180],'FontSize',24,'XTick',0:30:180)
legend({'GMR-Gal4/+','GMR-Gal4/ShiTS'},'Location','southwest')
legend('boxoff')
ylabel({'Normalized LDM'})
xlabel('Time (mins)')
shg
name = strcat(names(n),'_Sliding');
print(name{:},'-dpdf','-bestfit')


%% Generate Fig. 2A (Sliding of Shift)
n = [19];
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(any(logical(importedData.data(n,:)),1));
idx = ismember(upper(subArray(:,3)),ngn) & all(subData(:,[1 6]) > 0.03,2);
subArray2 = subArray(idx,:);


timebin = 80;
dEnd = [0:5:(180-timebin)]*60;
timebin = timebin*60;
timeMat = nan(size(subArray2,1),length(dEnd));
x = nan(size(subArray2,1), 2);
y = x;
nullX = x;
nullY = x;
ldmIso = nan(length(dEnd)-1,2); ldmShi = ldmIso;
proportion = nan(2,2,length(dEnd)-1);
for jj = 1:length(dEnd)-1
    for ii = 1:size(subArray2,1)
        turns = subArray2{ii,1};
        times = subArray2{ii,2};
        lights = subArray2{ii,11};
        
        idx = times >= dEnd(jj) & times < dEnd(jj) + timebin;
        if sum(idx) > turnfilter
        x(ii,:) = turnbias(turns(idx & lights));
        y(ii,:) = turnbias(turns(idx & ~lights));
        
        if jj == 1
            nullX(ii,:) = x(ii,:);
            nullY(ii,:) = y(ii,:);
        end
        end
    end
    ind= strcmp(subArray2(:,4),'ISO');

    dat = [nullX(ind,1),nullY(ind,1),x(ind,1), y(ind,1)];
    err = [nullX(ind,2),nullY(ind,2),x(ind,2), y(ind,2)];
    
    f = @(x) (x - nanmean(x,1))./nanstd(x,[],1);
    dataMat = f(dat);
    
    dataMat = dataMat(:,[3:4,1:2]);
    dataMat = dataMat(all(~isnan(dataMat),2),:);
   
    [~,temp] = min(abs(dataMat(:,1) - dataMat(:,3:4)),[],2);
    temp2 = temp == 1;
    
    [~,temp] =  min(abs(dataMat(:,2) - dataMat(:,3:4)),[],2);
    temp3 = temp == 1;
    
    temp4 = double([temp2,temp3]);
    
    temp4(temp4 == 0) = -1;
    
    proportion(1,:,jj) = nanmean(temp4,1);


    
    ind= strcmp(subArray2(:,4),'SHI');    
    dat = [nullX(ind,1),nullY(ind,1),x(ind,1), y(ind,1)];
    err = [nullX(ind,2),nullY(ind,2),x(ind,2), y(ind,2)];
    
    f = @(x) (x - nanmean(x,1))./nanstd(x,[],1);
    dataMat = f(dat);
    
    dataMat = dataMat(:,[3:4,1:2]);
   
    [~,temp] = min(abs(dataMat(:,1) - dataMat(:,3:4)),[],2);
    temp2 = temp == 1;
    
    [~,temp] =  min(abs(dataMat(:,2) - dataMat(:,3:4)),[],2);
    temp3 = temp == 1;
    
    temp4 = double([temp2,temp3]);
    
    temp4(temp4 == 0) = -1;
    
    proportion(2,:,jj) = nanmean(temp4,1);
    
    
end

clf
xVect = (dEnd(1:(end-1))+timebin/2)/60;
hPlot = plot(xVect,reshape(proportion,4,length(dEnd)-1)','LineWidth',3)

hPatch1 = patch([xVect fliplr(xVect)],[ldmIso(:,1)+ldmIso(:,2); flipud(ldmIso(:,1)-ldmIso(:,2))],...
    hPlot(1).Color);
hPatch2 = patch([xVect fliplr(xVect)],[ldmShi(:,1)+ldmShi(:,2); flipud(ldmShi(:,1)-ldmShi(:,2))],...
    hPlot(2).Color);

set([hPatch1,hPatch2],'FaceAlpha',0.25,'EdgeAlpha',0)
set(gca,'XLim',[0 180],'FontSize',24,'XTick',0:30:180)
legend({'GMR-Gal4/+','GMR-Gal4/ShiTS'},'Location','southwest')
legend('boxoff')
ylabel({'Normalized LDM'})
xlabel('Time (mins)')
shg
name = strcat(names(n),'_Sliding');
print(name{:},'-dpdf','-bestfit')

%% Make kymograph of individuals
n = [20];
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(any(logical(importedData.data(n,:)),1));
idx = ismember(upper(subArray(:,3)),ngn) & all(subData(:,[1 6]) > 0.03,2);


subArray2 = subArray(idx,:);
figure
ind= strcmp(subArray2(:,4),'ISO');

subArray2 = subArray2(ind,:);

timebin = 60;
dEnd = [0:10:(210-timebin)]*60;
timebin = timebin*60;
timeMat = nan(size(subArray2,1),length(dEnd));
x = nan(size(subArray2,1), 2);
y = x;
nullX = x;
nullY = x;
ldmIso = nan(length(dEnd)-1,2); ldmShi = ldmIso;
proportion = nan(2,2,length(dEnd)-1);
ldm = nan(size(subArray2,1),length(dEnd)-1)
for jj = 1:length(dEnd)
    for ii = 1:size(subArray2,1)
        turns = subArray2{ii,1};
        times = subArray2{ii,2};
        lights = subArray2{ii,11};
        
        idx = times >= dEnd(jj) & times < dEnd(jj) + timebin;
        if sum(idx) > turnfilter
        x = turnbias(turns(idx & lights));
        y = turnbias(turns(idx & ~lights));
        
        ldm(ii,jj) = x(:,1)-y(:,1);
        end
    end
end


dat = sortrows(abs(ldm),'descend')';
dat = dat(:,all(~isnan(dat),1));

xVect = (dEnd(1:(end-1)) + timebin/2)/60
subplot(2,1,2)
plot(xVect,smooth(nanmean(diff(dat,[],1),2)))
colormap(parula(1000))
caxis([-0.01 0.01])
pbaspect([1 1 1])
colorbar



subplot(2,1,1)
plot([dat(1,:);dat(end,:)]','.','MarkerSize',10)
pbaspect([1 1 1])
set(gca,'YLim',[0 0.5],'XLim',[0 200],'FontName','Arial')
xlabel('Individuals')
ylabel('LDM')
legend('Low Temp','High Temp')
shg

name = strcat(names(n),'_Individuals');
print(name{:},'-dpdf','-bestfit')


kde(dat(1,:))
hold on
kde(dat(end,:))

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

%% Make a sliding window of photoreceptor ldm over the whole experiment

neuron = ismember(upper(array(:,3)),'RH1');
idx1 = neuron & strcmp(effectors,'ISO');

idx2 = neuron & strcmp(effectors,'SHI');

%% Make individual plot

n = 7
    clf
%subplot(4,6,n)
behave = 2;
[ug,~,ugi] = uniqueRowsCA([upper(subArray(:,3))]);
ngn = genonames(logical(importedData.data(n,:)));
idx = ismember(upper(subArray(:,3)),ngn);
% split = 'SS02296'
%idx = ismember(array(:,3),split) 
effectors = subArray(:,4);
name = strcat(names{n},'_individuals')

nyidalur=interp1([1 128 129 256],[1 0 1; .2559 .248 .2559; .248 .2559 .248; 0 1 0],1:256);

[a,b] = sort(dat1 - dat3);

%[a,b] = sort(dat2 - dat4);

idx = b(end-2);

timeRange = [0 180]*60;
timeBins = [timeRange(1):900:timeRange(2)];
xVals = 7.5:15:240;

bias = nan(length(timeBins)-1,2);

tStamps = subArray{idx,2};
turns = subArray{idx,1};

turns = turns(tStamps>=timeRange(1) & tStamps<timeRange(2));
tStamps = tStamps(tStamps>=timeRange(1) & tStamps<timeRange(2));

turnFilter = 10;
for ii = 1:length(timeBins)-1
    tI = tStamps >= timeBins(ii) & tStamps < timeBins(ii+1);
    turn = turns(tI); 
    if length(turn) > turnFilter
    bias(ii,:) = turnbias(turn);
    end
end

figure
hPlot = plot(xVals,bias(:,1));
set(hPlot,'LineWidth',2,'LineStyle','--','Marker','s','MarkerSize',12,'Color','black',...
    'MarkerFaceColor','black')
set(hPlot.Parent,'YLim',[0 1])

hold on
xS = 10;
offset = 1.23;
for jj = 1:length(turns)
        x = [tStamps(jj)-xS,tStamps(jj)+xS,tStamps(jj)+xS,tStamps(jj)-xS]/60 - timeRange(1)/60;
    if turns(jj)
        y = [0 0 0.1 0.1] + offset;
        col = nyidalur(50,:);
    else
        y = [0 0 -0.1 -0.1] + offset;
        col = nyidalur(end-50,:);
    end
    hPatch = patch(x,y,col,'EdgeAlpha',0,'FaceAlpha',0.4);
end


lightV = [repmat([0.2 0.2 0 0],1,length(xVals)/2)];
y = interleave2(timeBins,timeBins)/60 - timeRange(1)/60;
offset3 = 1.5;
hPlot2 = plot(y(2:end-1),lightV+offset3,'Color','black','LineWidth',2);

yMax = 1.75;
x = interleave2(15:15:240,15:15:240)';
y = repmat([0 yMax yMax 0],1,floor(length(x)/4));
patch(x,y,[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0)

patch([0 250 250 0],[1 1 1.05 1.05],[1 1 1],'EdgeAlpha',0)
offset2 = 0.4;
patch([0 250 250 0],[1 1 1.05 1.05]+offset2,[1 1 1],'EdgeAlpha',0)



set(gca,'YLim',[0 yMax],'YTick',[0:0.5:1],'Box','off','FontSize',28,'FontWeight','bold',...
    'XLim',[0 250],'XTick',0:60:240,'TickDir','out')
ylabel('Turn Bias','FontSize',22)
xlabel('Time (mins)','FontSize',22)

set(gcf,'Renderer','painters')

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
filepath = '/Volumes/LaCie/Data/split_Rubin_Annotation.csv';
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
                VBE = bootstrp(1000,@(x,y) IBEV2(x,y),data(idx,[ii+10 ii+15]),...
                    error(idx,[ii+10 ii+15]));
                VBEMeanHigh{jj,kk}(ii) = nanmean(VBE);
                VBEErrorHigh{jj,kk}(ii) = nanstd(VBE);
                numFliesHigh(jj,kk) = sum(idx);
            end
            
            for ii = 1:5
                VBE = bootstrp(1000,@(x,y) IBEV2(x,y),data(idx,[ii ii+5]),...
                    error(idx,[ii+10 ii+15]));
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

%%
[a1,b1,c1] = uniqueRowsCA([subArray(:,3),subArray(:,4)]);

out = nan(length(a1),10)
for ii = 1:length(a1)
    idx = c1 == ii;
    subArray2 = subArray(idx,:);
    subData2 = subData(idx,:);
    [a2,b2,c2] = uniqueRowsCA([cellfun(@(x) x(2),subArray2(:,9)),...
        subArray2(:,7),cellfun(@(x) x(4),subArray2(:,9))]);
    temp = nan(size(a2,1),1);
    for kk = 1:size(a2,1)
    temp(kk,:) = mad(subData2(c2==kk,[2]) - subData2(c2==kk,[7]));
    end
    out(ii,1:length(temp)) = temp;
end
    
    

%subData

%% Make behavioral metric array
% make cell column name vector

cellColNames(13:42) = {'totalTurnsLow','lightTurnsLow','darkTurnsLow','totalTurnsHigh',...
    'lightTurnsHigh','darkTurnsHigh','totalTBLow','lightTBLow','darkTBLow','totalTBHigh',...
    'lightTBHigh','darkTBHigh','totalClumpLow','lightClumpLow','darkClumpLow','totalClumpHigh',...
    'lightClumpHigh','darkClumpHigh','totalSwitchLow','lightSwitchLow','darkSwitchLow',...
    'totalSwitchHigh','lightSwitchHigh','darkSwitchHigh','totalWallDistLow',...
    'lightWallDistLow','darkWallDistLow','totalWallDistHigh','lightWallDistHigh',...
    'darkWallDistHigh'};


cellColNames(43:60) = {'totalTBLowE','lightTBLowE','darkTBLowE','totalTBHighE',...
    'lightTBHighE','darkTBHighE','totalClumpLowE','lightClumpLowE','darkClumpLowE','totalClumpHighE',...
    'lightClumpHighE','darkClumpHighE','totalSwitchLowE','lightSwitchLowE','darkSwitchLowE',...
    'totalSwitchHighE','lightSwitchHighE','darkSwitchHighE'};

cellColNames(61:92) = {'ldmTurnsLow','ldmTurnsHigh','ldmTBLow','ldmTBHigh',...
    'ldmClumpLow','ldmClumpHigh','ldmSwitchLow','ldmSwitchHigh',...
    'ldmWallDistLow','ldmWallDistHigh','shiftTurnsLight','shiftTurnsDark',...
    'shiftTBLight','shiftTBDark','shiftClumpLight','shiftClumpDark',...
    'shiftSwitchLight','shiftSwitchDark','shiftWallDistLight','shiftWallDistDark',...
    'ldmTBLowE','ldmTBHighE',...
    'ldmClumpLowE','ldmClumpHighE','ldmSwitchLowE','ldmSwitchHighE',...
    'shiftTBLightE','shiftTBDarkE',...
    'shiftClumpLightE','shiftClumpDarkE','shiftSwitchLightE','shiftSwitchDarkE'};


cellColNames(93:182) = {'splitFirstTotalTurnsLow','splitSecondTotalTurnsLow',...
    'splitFirstLightTurnsLow','splitSecondLightTurnsLow',...
    'splitFirstDarkTurnsLow','splitSecondDarkTurnsLow',...
    'splitFirstTotalTurnsHigh','splitSecondTotalTurnsHigh',...
    'splitFirstLightTurnsHigh','splitSecondLightTurnsHigh',...
    'splitFirstDarkTurnsHigh','splitSecondDarkTurnsHigh',...
    'splitFirstTotalTBLow','splitSecondTotalTBLow',...
    'splitFirstLightTBLow','splitSecondLightTBLow',...
    'splitFirstDarkTBLow','splitSecondDarkTBLow',...
    'splitFirstTotalTBHigh','splitSecondTotalTBHigh',...
    'splitFirstLightTBHigh','splitSecondLightTBHigh',...
    'splitFirstDarkTBHigh','splitSecondDarkTBHigh',...
    'splitFirstTotalClumpLow','splitSecondTotalClumpLow',...
    'splitFirstLightClumpLow','splitSecondLightClumpLow',...
    'splitFirstDarkClumpLow','splitSecondDarkClumpLow',...
    'splitFirstTotalClumpHigh','splitSecondTotalClumpHigh',...
    'splitFirstLightClumpHigh','splitSecondLightClumpHigh',...
    'splitFirstDarkClumpHigh','splitSecondDarkClumpHigh',...
    'splitFirstTotalSwitchLow','splitSecondTotalSwitchLow',...
    'splitFirstLightSwitchLow','splitSecondLightSwitchLow',...
    'splitFirstDarkSwitchLow','splitSecondDarkSwitchLow',...
    'splitFirstTotalSwitchHigh','splitSecondTotalSwitchHigh',...
    'splitFirstLightSwitchHigh','splitSecondLightSwitchHigh',...
    'splitFirstDarkSwitchHigh','splitSecondDarkSwitchHigh',...
    'splitFirstTotalWallDistLow','splitSecondTotalWallDistLow',...
    'splitFirstLightWallDistLow','splitSecondLightWallDistLow',...
    'splitFirstDarkWallDistLow','splitSecondDarkWallDistLow',...
    'splitFirstTotalWallDistHigh','splitSecondTotalWallDistHigh',...
    'splitFirstLightWallDistHigh','splitSecondLightWallDistHigh',...
    'splitFirstDarkWallDistHigh','splitSecondDarkWallDistHigh',...
    'splitDiffTotalTurnsLow','splitDiffLightTurnsLow',...
    'splitDiffDarkTurnsLow','splitDiffTotalTurnsHigh',...
    'splitDiffLightTurnsHigh','splitDiffDarkTurnsHigh',...
    'splitDiffTotalTBLow','splitDiffLightTBLow',...
    'splitDiffDarkTBLow','splitDiffTotalTBHigh',...
    'splitDiffLightTBHigh','splitDiffDarkTBHigh',...
    'splitDiffTotalClumpLow','splitDiffLightClumpLow',...
    'splitDiffDarkClumpLow','splitDiffTotalClumpHigh',...
    'splitDiffLightClumpHigh','splitDiffDarkClumpHigh',...
    'splitDiffTotalSwitchLow','splitDiffLightSwitchLow',...
    'splitDiffDarkSwitchLow','splitDiffTotalSwitchHigh',...
    'splitDiffLightSwitchHigh','splitDiffDarkSwitchHigh',...
    'splitDiffTotalWallDistLow','splitDiffLightWallDistLow',...
    'splitDiffDarkWallDistLow','splitDiffTotalWallDistHigh',...
    'splitDiffLightWallDistHigh','splitDiffDarkWallDistHigh'};




% time in minutes
lowT = [0 60]*60;
highT = [120 200]*60;

% threshold
turnThreshold = 100;
numRe = 1;

% other parameters
armWidth = 3; % in mm

% turnBias calculation
for ii = 1:size(allDataCell,1)
    
    if mod(ii,100) == 0
        ii
    end
    
    % return index for all low and high temp turns
    idxL = allDataCell{ii,2} >= lowT(1) & allDataCell{ii,2} < lowT(2);
    idxH = allDataCell{ii,2} >= highT(1) & allDataCell{ii,2} < highT(2);
    
    %idxH(cumsum(idxH) > 200) = 0;
    
    % skip the fly if below threshold
    if sum(idxH) < turnThreshold || sum(idxL) < turnThreshold
        continue
    end
    
    % total turns low
    allDataCell(ii,13) = num2cell(length(allDataCell{ii,1}(idxL)));
    % light turns low
    allDataCell(ii,14) = num2cell(length(allDataCell{ii,1}(idxL & allDataCell{ii,11})));
    % dark turns low
    allDataCell(ii,15) = num2cell(length(allDataCell{ii,1}(idxL & ~allDataCell{ii,11})));
    
    
    % total turns high
    allDataCell(ii,16) = num2cell(length(allDataCell{ii,1}(idxH)));
    % light turns high
    allDataCell(ii,17) = num2cell(length(allDataCell{ii,1}(idxH & allDataCell{ii,11})));
    % dark turns high
    allDataCell(ii,18) = num2cell(length(allDataCell{ii,1}(idxH & ~allDataCell{ii,11})));
    
    % total tb low
    allDataCell(ii,19) = num2cell(nanmean(allDataCell{ii,1}(idxL)));
    % light tb low
    allDataCell(ii,20) = num2cell(nanmean(allDataCell{ii,1}(idxL & allDataCell{ii,11})));
    % dark tb low
    allDataCell(ii,21) = num2cell(nanmean(allDataCell{ii,1}(idxL & ~allDataCell{ii,11})));
    
    % total tb high
    allDataCell(ii,22) = num2cell(nanmean(allDataCell{ii,1}(idxH)));
    % light tb high
    allDataCell(ii,23) = num2cell(nanmean(allDataCell{ii,1}(idxH & allDataCell{ii,11})));
    % dark tb high
    allDataCell(ii,24) = num2cell(nanmean(allDataCell{ii,1}(idxH & ~allDataCell{ii,11})));
    
    % clumpiness
    
    % total clump low
    time = lowT(2) - lowT(1);
    tstamps = allDataCell{ii,2}(idxL);
    turns = allDataCell{ii,1}(idxL);
    iti = diff(allDataCell{ii,2}(idxL));
    allDataCell(ii,25) = num2cell(mad(iti,1)/(time/length(turns)));
    % light clump low
    tstamps = allDataCell{ii,2}(idxL & ~allDataCell{ii,11});
    turns = allDataCell{ii,1}(idxL & allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxL & allDataCell{ii,11}));
    allDataCell(ii,26) = num2cell(mad(iti,1)/(time/length(turns)));
    % dark clump low
    turns = allDataCell{ii,1}(idxL & ~allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxL & ~allDataCell{ii,11}));
    allDataCell(ii,27) = num2cell(mad(iti,1)/(time/length(turns)));
    
    % total clump high
    time = highT(2) - highT(1);
    turns = allDataCell{ii,1}(idxH);
    iti = diff(allDataCell{ii,2}(idxH));
    allDataCell(ii,28) = num2cell(mad(iti,1)/(time/length(turns)));
    % light clump high
    turns = allDataCell{ii,1}(idxH & allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxH & allDataCell{ii,11}));
    allDataCell(ii,29) = num2cell(mad(iti,1)/(time/length(turns)));
    % dark clump high
    turns = allDataCell{ii,1}(idxH & ~allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxH & ~allDataCell{ii,11}));
    allDataCell(ii,30) = num2cell(mad(iti,1)/(time/length(turns)));
    
    % switchiness
    
    % total switch low
    turns = allDataCell{ii,1}(idxL);
    iti = diff(allDataCell{ii,2}(idxL));
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,31} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    allDataCell{ii,183} = nanmean(iti(turnLogical));
    allDataCell{ii,184} = nanmean(iti(~turnLogical));
    % light switch low
    turns = allDataCell{ii,1}(idxL & allDataCell{ii,11});
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,32} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    % dark switch low
    turns = allDataCell{ii,1}(idxL & ~allDataCell{ii,11});
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,33} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    
    % total switch high
    turns = allDataCell{ii,1}(idxH);
    iti = diff(allDataCell{ii,2}(idxH));
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,34} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    allDataCell{ii,185} = nanmean(iti(turnLogical));
    allDataCell{ii,186} = nanmean(iti(~turnLogical));
    % light switch high
    turns = allDataCell{ii,1}(idxH & allDataCell{ii,11});
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,35} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    % dark switch high
    turns = allDataCell{ii,1}(idxH & ~allDataCell{ii,11});
    turnLogical = (turns(1:(end-1))==1 & turns(2:end))==1;
    allDataCell{ii,36} = sum(turnLogical)/(sum(turns)*nanmean(turns));
    
    % wall following low
    
    turns = allDataCell{ii,1};
    % total wall low left
    wallsL = allDataCell{ii,10}(idxL & ~turns);
    wallsR = allDataCell{ii,10}(idxL & turns);
    allDataCell{ii,37} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    % light wall low
    wallsL = allDataCell{ii,10}(idxL & ~turns & allDataCell{ii,11});
    wallsR = allDataCell{ii,10}(idxL & turns & allDataCell{ii,11});
    allDataCell{ii,38} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    % dark wall low left
    wallsL = allDataCell{ii,10}(idxL & ~turns & ~allDataCell{ii,11});
    wallsR = allDataCell{ii,10}(idxL & turns & ~allDataCell{ii,11});
    allDataCell{ii,39} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    
    % wall following high
    
    turns = allDataCell{ii,1};
    % total wall low left
    walls = allDataCell{ii,10}(idxH)%; & ~turns);
    wallsR = allDataCell{ii,10}(idxH & turns);
    allDataCell{ii,40} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    % light wall low
    wallsL = allDataCell{ii,10}(idxH & ~turns & allDataCell{ii,11});
    wallsR = allDataCell{ii,10}(idxH & turns & allDataCell{ii,11});
    allDataCell{ii,41} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    % dark wall low left
    wallsL = allDataCell{ii,10}(idxH & ~turns & ~allDataCell{ii,11});
    wallsR = allDataCell{ii,10}(idxH & turns & ~allDataCell{ii,11});
    allDataCell{ii,42} = (1-nanmean([(1-wallsL);wallsR]))*armWidth;
    
    
    
    % calculate the errors for each of the previous stats
    
    % error on tb
    
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'totalTBLow')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'totalTurnsLow')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,43} = error;
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBLow')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsLow')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,44} = error;
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBLow')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsLow')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,45} = error;
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'totalTBHigh')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'totalTurnsHigh')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,46} = error;
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBHigh')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsHigh')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,47} = error;
    % total tb
    p = cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBHigh')));
    n = cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsHigh')));
    error = sqrt(n.*p.*(1-p))./n;
    allDataCell{ii,48} = error;
    
    % return index for all low and high temp turns
    idxL = allDataCell{ii,2} >= lowT(1) & allDataCell{ii,2} < lowT(2);
    idxH = allDataCell{ii,2} >= highT(1) & allDataCell{ii,2} < highT(2);
    
    
    % error on clump
    time = lowT(2) - lowT(1);
    % total low temp
    turns = allDataCell{ii,1}(idxL);
    iti = diff(allDataCell{ii,2}(idxL));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,49} = nanstd(clumpRe);
    % light low temp
    turns = allDataCell{ii,1}(idxL & allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxL & allDataCell{ii,11}));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,50} = nanstd(clumpRe);
    % dark low temp
    turns = allDataCell{ii,1}(idxL & ~allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxL & ~allDataCell{ii,11}));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,51} = nanstd(clumpRe);
    % total high temp
    turns = allDataCell{ii,1}(idxH);
    iti = diff(allDataCell{ii,2}(idxH));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,52} = nanstd(clumpRe);
    % light high temp
    turns = allDataCell{ii,1}(idxH & allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxH & allDataCell{ii,11}));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,53} = nanstd(clumpRe);
    % dark high temp
    turns = allDataCell{ii,1}(idxH & ~allDataCell{ii,11});
    iti = diff(allDataCell{ii,2}(idxH & ~allDataCell{ii,11}));
    clumpRe = nan(100,1);
    for jj = 1:numRe
        itiRe = randsample(iti,length(turns),true);
        clumpRe(jj) = mad(itiRe,1)/(time/length(turns));
    end
    allDataCell{ii,54} = nanstd(clumpRe);
    
    % error on switch
    
    % total low
    turns = allDataCell{ii,1}(idxL);
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,55} = nanstd(swiRe);
    % light low
    turns = allDataCell{ii,1}(idxL & allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,56} = nanstd(swiRe);
    % dark low
    turns = allDataCell{ii,1}(idxL & ~allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,57} = nanstd(swiRe);
    % total high
    turns = allDataCell{ii,1}(idxH);
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,58} = nanstd(swiRe);
    % light high
    turns = allDataCell{ii,1}(idxH & allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,59} = nanstd(swiRe);
    % dark high
    turns = allDataCell{ii,1}(idxH & ~allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    swi = turns(1:end-1)+turns(2:end);
    for jj = 1:numRe
        swiIdx = randi(length(swi),length(swi),1);
        turnRe = turns(swiIdx);
        divRe = 2*length(turnRe)*nanmean(turnRe)*(1-nanmean(turnRe));
        swiRe(jj) = sum(swi(swiIdx) == 1)/div;
    end
    allDataCell{ii,60} = nanstd(swiRe);
    
    % ldm turns for each stat
    % ldm turns low
    allDataCell(ii,61) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsLow')))};
    % ldm turns high
    allDataCell(ii,62) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsHigh'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsHigh')))};
    
    % ldm tb for each stat
    % ldm tb low
    allDataCell(ii,63) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBLow')))};
    % ldm tb high
    allDataCell(ii,64) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBHigh'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBHigh')))};
    
    % ldm clump for each stat
    % ldm low
    allDataCell(ii,65) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpLow')))};
    % ldm high
    allDataCell(ii,66) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpHigh'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpHigh')))};
    
    % ldm switch for each stat
    % ldm low
    allDataCell(ii,67) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchLow')))};
    % ldm tb high
    allDataCell(ii,68) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchHigh'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchHigh')))};
    
    % ldm wall dist for each stat
    % ldm low
    allDataCell(ii,69) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightWallDistLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkWallDistLow')))};
    % ldm tb high
    allDataCell(ii,70) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightWallDistHigh'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkWallDistHigh')))};
    
    
    % shift low to high light
    allDataCell(ii,71) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTurnsHigh')))};
    % shift low to high dark
    allDataCell(ii,72) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTurnsHigh')))};
    
    
    % shift low to high light
    allDataCell(ii,73) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBHigh')))};
    % shift low to high dark
    allDataCell(ii,74) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBHigh')))};
    
    
    % shift low to high light
    allDataCell(ii,75) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpHigh')))};
    % shift low to high dark
    allDataCell(ii,76) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpHigh')))};
    
    
    % shift low to high light
    allDataCell(ii,77) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchHigh')))};
    % shift low to high dark
    allDataCell(ii,78) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchHigh')))};
    
    
    % shift low to high light
    allDataCell(ii,79) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'lightWallDistLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightWallDistHigh')))};
    % shift low to high dark
    allDataCell(ii,80) = {cell2mat(allDataCell(ii,strcmp(cellColNames,'darkWallDistLow'))) - ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkWallDistHigh')))};
    
    % ldm tb error for each stat
    % ldm low
    allDataCell(ii,81) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBLowE'))).^2)};
    % ldm tb high
    allDataCell(ii,82) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBHighE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBHighE'))).^2)};
    % ldm clump error for each stat
    % ldm low
    allDataCell(ii,83) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpLowE'))).^2)};
    % ldm tb high
    allDataCell(ii,84) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpHighE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpHighE'))).^2)};
    % ldm tb error for each stat
    % ldm low
    allDataCell(ii,85) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchLowE'))).^2)};
    % ldm tb high
    allDataCell(ii,86) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchHighE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchHighE'))).^2)};
    
    % shift low to high light
    allDataCell(ii,87) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightTBHighE'))).^2)};
    % shift low to high dark
    allDataCell(ii,88) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkTBHighE'))).^2)};
    
    
    % shift low to high light
    allDataCell(ii,89) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightClumpHighE'))).^2)};
    % shift low to high dark
    allDataCell(ii,90) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkClumpHighE'))).^2)};
    
    
    % shift low to high light
    allDataCell(ii,91) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'lightSwitchHighE'))).^2)};
    % shift low to high dark
    allDataCell(ii,92) = {sqrt(cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchLowE'))).^2 + ...
        cell2mat(allDataCell(ii,strcmp(cellColNames,'darkSwitchHighE'))).^2)};
    
    
    
    % return index for all low and high temp turns first half and second half of each time block
    idxL1 = allDataCell{ii,2} >= lowT(1) & allDataCell{ii,2} < (lowT(2)-lowT(1))/2+lowT(1);
    idxL2 = allDataCell{ii,2} >= (lowT(2)-lowT(1))/2+lowT(1) & allDataCell{ii,2} < lowT(2);
    idxH1 = allDataCell{ii,2} >= highT(1) & allDataCell{ii,2} < (highT(2)-highT(1))/2+highT(1);
    idxH2 = allDataCell{ii,2} >= (highT(2)-highT(1))/2+highT(1) & allDataCell{ii,2} < highT(2);
    
    
    % split turns low
    allDataCell{ii,93} = length(allDataCell{ii,1}(idxL1));
    allDataCell{ii,94} = length(allDataCell{ii,1}(idxL2));
    allDataCell{ii,153} = allDataCell{ii,94} - allDataCell{ii,93};
    % light turns low
    allDataCell{ii,95} = length(allDataCell{ii,1}(idxL1 & allDataCell{ii,11}));
    allDataCell{ii,96} = length(allDataCell{ii,1}(idxL2 & allDataCell{ii,11}));
    allDataCell{ii,154} = allDataCell{ii,96} - allDataCell{ii,95};
    % dark turns low
    allDataCell{ii,97} = length(allDataCell{ii,1}(idxL1 & allDataCell{ii,11}));
    allDataCell{ii,98} = length(allDataCell{ii,1}(idxL2 & allDataCell{ii,11}));
    allDataCell{ii,155} = allDataCell{ii,98} - allDataCell{ii,97};
    
    % split turns high
    allDataCell{ii,99} = length(allDataCell{ii,1}(idxH1));
    allDataCell{ii,100} = length(allDataCell{ii,1}(idxH2));
    allDataCell{ii,156} = allDataCell{ii,100} - allDataCell{ii,99};
    % light turns high
    allDataCell{ii,101} = length(allDataCell{ii,1}(idxH1 & allDataCell{ii,11}));
    allDataCell{ii,102} = length(allDataCell{ii,1}(idxH2 & allDataCell{ii,11}));
    allDataCell{ii,157} = allDataCell{ii,102} - allDataCell{ii,101};
    % dark turns high
    allDataCell{ii,103} = length(allDataCell{ii,1}(idxH1 & allDataCell{ii,11}));
    allDataCell{ii,104} = length(allDataCell{ii,1}(idxH2 & allDataCell{ii,11}));
    allDataCell{ii,158} = allDataCell{ii,104} - allDataCell{ii,103};
    
    % total tb low
    allDataCell{ii,105} = nanmean(allDataCell{ii,1}(idxL1));
    allDataCell{ii,106} = nanmean(allDataCell{ii,1}(idxL2));
    allDataCell{ii,159} = allDataCell{ii,105} - allDataCell{ii,106};
    % light tb low
    allDataCell{ii,107} = nanmean(allDataCell{ii,1}(idxL1 & allDataCell{ii,11}));
    allDataCell{ii,108} = nanmean(allDataCell{ii,1}(idxL2 & allDataCell{ii,11}));
    allDataCell{ii,160} = allDataCell{ii,108} - allDataCell{ii,107};
    % dark tb low
    allDataCell{ii,109} = nanmean(allDataCell{ii,1}(idxL1 & ~allDataCell{ii,11}));
    allDataCell{ii,110} = nanmean(allDataCell{ii,1}(idxL2 & ~allDataCell{ii,11}));
    allDataCell{ii,161} = allDataCell{ii,110} - allDataCell{ii,109};
    
    % total tb high
    allDataCell{ii,111} = nanmean(allDataCell{ii,1}(idxH1));
    allDataCell{ii,112} = nanmean(allDataCell{ii,1}(idxH2));
    allDataCell{ii,162} = allDataCell{ii,112} - allDataCell{ii,111};
    % light tb high
    allDataCell{ii,113} = nanmean(allDataCell{ii,1}(idxH1 & allDataCell{ii,11}));
    allDataCell{ii,114} = nanmean(allDataCell{ii,1}(idxH2 & allDataCell{ii,11}));
    allDataCell{ii,163} = allDataCell{ii,114} - allDataCell{ii,113};
    % dark tb high
    allDataCell{ii,115} = nanmean(allDataCell{ii,1}(idxH1 & ~allDataCell{ii,11}));
    allDataCell{ii,116} = nanmean(allDataCell{ii,1}(idxH2 & ~allDataCell{ii,11}));
    allDataCell{ii,164} = allDataCell{ii,116} - allDataCell{ii,115};
    
    
    % split clumpiness
    
    % total low
    time = (lowT(2)-lowT(1))/2;
    turns1 = allDataCell{ii,1}(idxL1);
    turns2 = allDataCell{ii,1}(idxL2);
    iti1 = diff(allDataCell{ii,2}(idxL1));
    iti2 = diff(allDataCell{ii,2}(idxL2));
    allDataCell{ii,117} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,118} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,165} = allDataCell{ii,118} - allDataCell{ii,117};
    
    % light clump low
    turns1 = allDataCell{ii,1}(idxL1 & allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxL2 & allDataCell{ii,11});
    iti1 = diff(allDataCell{ii,2}(idxL1));
    iti2 = diff(allDataCell{ii,2}(idxL2));
    allDataCell{ii,119} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,120} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,166} = allDataCell{ii,120} - allDataCell{ii,119};
    % dark clump low
    turns1 = allDataCell{ii,1}(idxL1 & ~allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxL2 & ~allDataCell{ii,11});
    iti1 = diff(allDataCell{ii,2}(idxL1));
    iti2 = diff(allDataCell{ii,2}(idxL2));
    allDataCell{ii,121} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,122} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,167} = allDataCell{ii,122} - allDataCell{ii,121};
    % total high
    time = (highT(2)-highT(1))/2;
    turns1 = allDataCell{ii,1}(idxH1);
    turns2 = allDataCell{ii,1}(idxH2);
    iti1 = diff(allDataCell{ii,2}(idxH1));
    iti2 = diff(allDataCell{ii,2}(idxH2));
    allDataCell{ii,123} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,124} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,168} = allDataCell{ii,124} - allDataCell{ii,123};
    % light high
    turns1 = allDataCell{ii,1}(idxH1 & allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxH2 & allDataCell{ii,11});
    iti1 = diff(allDataCell{ii,2}(idxH1));
    iti2 = diff(allDataCell{ii,2}(idxH2));
    allDataCell{ii,125} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,126} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,169} = allDataCell{ii,126} - allDataCell{ii,125};
    % dark high
    turns1 = allDataCell{ii,1}(idxH1 & ~allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxH2 & ~allDataCell{ii,11});
    iti1 = diff(allDataCell{ii,2}(idxH1));
    iti2 = diff(allDataCell{ii,2}(idxH2));
    allDataCell{ii,127} = mad(iti1,1)/(time/length(turns1));
    allDataCell{ii,128} = mad(iti2,1)/(time/length(turns2));
    allDataCell{ii,170} = allDataCell{ii,128} - allDataCell{ii,127};
    
    % switchiness
    
    % total switch low
    turns1 = allDataCell{ii,1}(idxL1);
    turns2 = allDataCell{ii,1}(idxL2);
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,129} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,130} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,171} = allDataCell{ii,130} - allDataCell{ii,129};
    % light switch low
    turns1 = allDataCell{ii,1}(idxL1 & allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxL2 & allDataCell{ii,11});
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,131} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,132} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,172} = allDataCell{ii,132} - allDataCell{ii,131};
    % dark switch low
    turns1 = allDataCell{ii,1}(idxL1 & ~allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxL2 & ~allDataCell{ii,11});
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,133} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,134} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,173} = allDataCell{ii,134} - allDataCell{ii,133};
    
    % total switch high
    turns1 = allDataCell{ii,1}(idxH1);
    turns2 = allDataCell{ii,1}(idxH2);
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,135} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,136} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,174} = allDataCell{ii,136} - allDataCell{ii,135};
    % light switch high
    turns1 = allDataCell{ii,1}(idxH1 & allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxH2 & allDataCell{ii,11});
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,137} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,138} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,175} = allDataCell{ii,138} - allDataCell{ii,137};
    % dark switch high
    turns1 = allDataCell{ii,1}(idxH1 & ~allDataCell{ii,11});
    turns2 = allDataCell{ii,1}(idxH2 & ~allDataCell{ii,11});
    div1 = 2*length(turns1)*nanmean(turns1)*(1-nanmean(turns1));
    div2 = 2*length(turns2)*nanmean(turns2)*(1-nanmean(turns2));
    allDataCell{ii,139} = sum((turns1(1:end-1)+turns1(2:end))==1)/div1;
    allDataCell{ii,140} = sum((turns2(1:end-1)+turns2(2:end))==1)/div2;
    allDataCell{ii,176} = allDataCell{ii,140} - allDataCell{ii,139};
    
    % wall following low
    
    turns = allDataCell{ii,1};
    % total wall low
    wallsL1 = allDataCell{ii,10}(idxL1 & ~turns);
    wallsL2 = allDataCell{ii,10}(idxL2 & ~turns);
    wallsR1 = allDataCell{ii,10}(idxL1 & turns);
    wallsR2 = allDataCell{ii,10}(idxL2 & turns);
    allDataCell{ii,141} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,142} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,177} = allDataCell{ii,142} - allDataCell{ii,141};
    % light wall low
    wallsL1 = allDataCell{ii,10}(idxL1 & ~turns & allDataCell{ii,11});
    wallsL2 = allDataCell{ii,10}(idxL2 & ~turns & allDataCell{ii,11});
    wallsR1 = allDataCell{ii,10}(idxL1 & turns & allDataCell{ii,11});
    wallsR2 = allDataCell{ii,10}(idxL2 & turns & allDataCell{ii,11});
    allDataCell{ii,143} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,144} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,178} = allDataCell{ii,144} - allDataCell{ii,143};
    % dark wall low
    wallsL1 = allDataCell{ii,10}(idxL1 & ~turns & ~allDataCell{ii,11});
    wallsL2 = allDataCell{ii,10}(idxL2 & ~turns & ~allDataCell{ii,11});
    wallsR1 = allDataCell{ii,10}(idxL1 & turns & ~allDataCell{ii,11});
    allDataCell{ii,145} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,146} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,179} = allDataCell{ii,146} - allDataCell{ii,145};
    
    % wall following high
    
    turns = allDataCell{ii,1};
    % total wall low
    wallsL1 = allDataCell{ii,10}(idxH1 & ~turns);
    wallsL2 = allDataCell{ii,10}(idxH2 & ~turns);
    wallsR1 = allDataCell{ii,10}(idxH1 & turns);
    wallsR2 = allDataCell{ii,10}(idxH2 & turns);
    allDataCell{ii,147} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,148} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,180} = allDataCell{ii,148} - allDataCell{ii,147};
    % light wall low
    wallsL1 = allDataCell{ii,10}(idxH1 & ~turns & allDataCell{ii,11});
    wallsL2 = allDataCell{ii,10}(idxH2 & ~turns & allDataCell{ii,11});
    wallsR1 = allDataCell{ii,10}(idxH1 & turns & allDataCell{ii,11});
    wallsR2 = allDataCell{ii,10}(idxH2 & turns & allDataCell{ii,11});
    allDataCell{ii,149} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,150} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,181} = allDataCell{ii,150} - allDataCell{ii,149};
    % dark wall low
    wallsL1 = allDataCell{ii,10}(idxH1 & ~turns & ~allDataCell{ii,11});
    wallsL2 = allDataCell{ii,10}(idxH2 & ~turns & ~allDataCell{ii,11});
    wallsR1 = allDataCell{ii,10}(idxH1 & turns & ~allDataCell{ii,11});
    wallsR2 = allDataCell{ii,10}(idxH2 & turns & ~allDataCell{ii,11});
    allDataCell{ii,151} = (1-nanmean([(1-wallsL1);wallsR1]))*armWidth;
    allDataCell{ii,152} = (1-nanmean([(1-wallsL2);wallsR2]))*armWidth;
    allDataCell{ii,182} = allDataCell{ii,152} - allDataCell{ii,151};
    
    
end

% take only non-empty cells
allDataCell = allDataCell(~sum(cellfun(@isempty,allDataCell),2) > 0,:);


save('allBehaveMat.mat','allDataCell')
save('allBehaveColNames.mat','cellColNames')

%% If loading behavioral data

load('allBehaveMat.mat');

%% Calculate aggregate statistics

behaveStatsIdx = [13:42];
behaveStatsIdx(1:3:end) = [];
behaveStatsIdx = [20:21 63];
behaveMat = cell2mat(allDataCell(string(cell2mat(allDataCell(:,4))) == 'ISO',...
    behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1,'FontSize',16)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(jet)
colorbar
caxis([-1 1]);

savefig('ldmBehaviors')

%% Make genotype figure

behaveIdx = [13:42 61:70];
%behaveIdx = [20:21 23:24 63:64];
% behaveIdx = behaveIdx([...
%     1 2 5 6 9 10 13 14 17 18,3 4 7 8 11 12 15 16 19 20]);
behaveMat = cell2mat(allDataCell(:,behaveIdx,:));
[a1,b1,c1] = unique(allDataCell(:,3));
[a2,b2,c2] = unique(string(cell2mat(allDataCell(:,4)))); %unique(allDataCell(:,4));

behaveErrIdx = [43:60];
behaveErrMat = cell2mat(allDataCell(:,behaveErrIdx));

names = cellColNames(behaveIdx);
clear highs
for ii = 1:length(behaveIdx)
    highs(ii) = contains(names(ii),'High');
end
clear data error
for ii = 1:length(a1)
    for jj = 1:length(a2)
        idx = c1 == ii & c2 == jj;
        data(ii,:,jj) = nanmean(behaveMat(idx,:));
        error(ii,[1:6 25:40],jj) = nanstd(behaveMat(idx,[1:6 25:40]));
        error(ii,7:24,jj) = sqrt(abs(nanvar(behaveMat(idx,7:24))-nanmean(behaveErrMat(idx,:).^2)));
        flycount(ii,:,jj) = sum(idx);
        dataCorr{ii,:,jj} = behaveMat(idx,:);
    end
    
end

for kk = 1:numRe
    for ii = 1:length(a1)
        for jj = 1:length(a2)
            idx = find(c1 == ii & c2 == jj);
            idx = randsample(idx,length(idx),true);
            dataRe(ii,:,jj,kk) = nanmean(behaveMat(idx,:));
            errorRe(ii,[1:6 25:40],jj,kk) = nanstd(behaveMat(idx,[1:6 25:40]));
            errorRe(ii,7:24,jj,kk) = sqrt(abs(nanvar(behaveMat(idx,7:24))-nanmean(behaveErrMat(idx,:).^2)));
            
        end
    end
end

compIdx = 2;

data2 = data(:,highs,:)./data(:,~highs,:);
deltaData = (data2(:,:,compIdx) - data2(:,:,1))./data2(:,:,1);
finalIdx = sum(~isnan(deltaData),2) > 0;
deltaData = deltaData(finalIdx,:);

dataSE = nanstd(dataRe,[],4);
dataSE2 = dataSE(:,highs,:)./dataSE(:,~highs,:);
deltaDataSE = sqrt(dataSE2(:,:,compIdx).^2 + dataSE2(:,:,1).^2);
deltaDataSE = deltaDataSE(finalIdx,:);

cv = error;
cv2 = cv(:,highs,:)./cv(:,~highs,:);
deltaCV = (cv2(:,:,compIdx) - cv2(:,:,1))./cv2(:,:,1);
deltaCV = deltaCV(finalIdx,:);

CVSE = nanstd(errorRe,[],4);
CVSE2 = CVSE(:,highs,:)./CVSE(:,~highs,:);
deltaCVSE = sqrt(CVSE2(:,:,compIdx).^2 + CVSE2(:,:,1).^2);
deltaCVSE = deltaCVSE(finalIdx,:);

for ii = 1:length(dataCorr)
    if ~isempty(dataCorr{ii,:,compIdx})
        co = corr([dataCorr{ii,:,compIdx}(:,highs) dataCorr{ii,:,compIdx}(:,~highs)],...
            'rows','pairwise')
        dataCorrs(ii,:) = diag(co,20);
    end
end

for ii = 1:length(dataCorr)
    if ~isempty(dataCorr{ii,:,1})
        co = corr([dataCorr{ii,:,1}(:,highs) dataCorr{ii,:,1}(:,~highs)],...
            'rows','pairwise')
        dataCorrsIso(ii,:) = diag(co,20);
    end
end

gNames = a1(finalIdx);

% figure
% imagesc(deltaData)

finalNames = names(highs);
% for ii = 1:length(finalNames)
%     finalNames{ii} = finalNames{ii}(1:(end-4));
% end
% set(gca,'YTickLabel',gNames,'XTickLabel',finalNames,'XTickLabelRotation',45,...
%     'FontWeight','bold','YTick',1:length(a1),'XTick',1:length(behaveIdx))
% colormap(skalafell)
% colorbar
% caxis([-1 1]);
% shg
%
% % figure
% % imagesc(deltaCV)
%
% finalNames = names(highs);
% for ii = 1:length(finalNames)
%     finalNames{ii} = finalNames{ii}(1:(end-4));
% end
% set(gca,'YTickLabel',a1,'XTickLabel',finalNames,'XTickLabelRotation',45,...
%     'FontWeight','bold','YTick',1:length(a1),'XTick',1:length(behaveIdx))
% colormap(skalafell)
% colorbar
% caxis([-1 1]);
% shg
%
%%
cvCols = [5 6 17];
dataCols = [2 3 8 9 11 12 14 15];

data = (data(:,highs,:)-data(:,~highs,:))./data(:,~highs,:);
cv = (cv(:,highs,:)-cv(:,~highs,:))./cv(:,~highs,:);

effIdx = 2;
idx = ~sum(isnan(data(:,:,effIdx)),2) > 0;
data = data(idx,:,effIdx);
cv = cv(idx,:,effIdx);

plotData = [cv(:,cvCols) data(:,dataCols)];

gNames = a1(idx);

%plotData = [deltaCV(:,[5 6 17]) deltaData(:,dataCols)]

finalNames = names(highs);
for ii = 1:length(finalNames)
    finalNames{ii} = finalNames{ii}(1:(end-4));
end

idx = contains(gNames,'SS');

gNames = gNames(idx);

x = nanmean(plotData(idx,:));
y = nanstd(plotData(idx,:));
zscores = (plotData(idx,:) - x)./y;


corrMat = corr(zscores','rows','pairwise');
%%

metric = 'correlation'
Y = pdist(plotData(idx,:),metric);
Y = squareform(Y);
method = 'average'
Z = linkage(Y,method,metric);
[~,~,i] = dendrogram(Z,size(Z,1)+1);
%[~,i] = sort(plotData(:,2));
b = i;


genotypes = gNames(b);

corrMat = corrMat(i,i);
%subplot(1,3,3)
imagesc(corrMat)

set(gca,'YTickLabel',genotypes,'FontWeight','bold','YTick',1:length(genotypes),...
    'XTickLabel',genotypes,...
    'XTickLabelRotation',45,'XTick',1:length(genotypes),'FontSize',4)
colormap(gcf,jet)
colorbar
caxis([-1 1]);
pbaspect([1 1 1])


subplot(1,2,1)
imagesc(zscores(i,:))

set(gca,'YTickLabel',genotypes,'FontWeight','bold','YTick',1:length(genotypes),...
    'XTickLabel',finalNames([cvCols dataCols]),...
    'XTickLabelRotation',45,'XTick',1:(length(finalNames)*2),'FontSize',10)
colormap(gcf,skalafell)
colorbar
caxis([-3 3])
shg

filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

names = importedData.textdata(2:end,1);

genonames = importedData.textdata(1,2:end);

shorterNeuropilNames = {'PB','FB','EB','NO','Rubus','Crepine','Gall','LAL','PS'};

shorterNames = {'P-FN1/2','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','PF-FGall','PF-FRub',...
    'PF-LCre','LPs-P',...
    'PBCap','PBInt','Ps-P1','Ps-P2'};

rNamesSub = genotypes;

clear neuronLogicalMat
for ii = 1:length(genotypes)
    idx = strcmp(rNamesSub(ii),genonames);
    if sum(idx) > 0
        neuronLogicalMat(ii,:) = importedData.data(:,idx)';
        
    end
end

numClust = sum(sum(neuronLogicalMat(1:(end-1),:) == 1 & neuronLogicalMat(2:(end),:) == 1));
oldNumClust = nan(1);
bestMethod = [metric '_' method];
if numClust > oldNumClust
    bestMethod = [metric '_' method];
    oldNumClust = numClust;
end
bestMethod
subplot(1,2,2)
imagesc(neuronLogicalMat)

set(gca,'YTickLabel',rNamesSub,'FontWeight','bold','YTick',1:length(b),...
    'XTickLabel',names,'XTickLabelRotation',45,'XTick',1:length(behaveIdx),...
    'FontSize',10)
shg


%% Distribution of average iti's depending on right->right or right->left turns
rtr = cell2mat(allDataCell(:,185))
rtl = cell2mat(allDataCell(:,186))
clf
histogram(rtr,'Normalization','probability','BinEdges',0:0.1:30)
hold on
histogram(rtl,'Normalization','probability','BinEdges',0:0.1:30,'FaceAlpha',0.5)

% scatter(rtr,rtl)
% ylim([0 50])
% xlim([0 50])

%% Make some quick summary figures
behaveStatsIdx = [15:3:42];
%behaveStatsIdx(1:3:end) = []
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('darkBehavior')
%%
behaveStatsIdx = [13:42];
behaveStatsIdx(1:3:end) = [];
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('allBehaviorLightDark')

%%
behaveStatsIdx = [13:42];
behaveStatsIdx(1:3:end) = [];
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('allBehaviorLightDark')

%%
behaveStatsIdx = [93:152];
behaveStatsIdx(1:6:end) = [];
behaveStatsIdx(1:5:end) = [];
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('allBehaviorLightDarkIntraindividual')

%%
behaveStatsIdx = [61:70];
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('ldmBehaviors')


%%
behaveStatsIdx = [71:80];
%61:80 93:182]
behaveMat = cell2mat(allDataCell(:,behaveStatsIdx));
x = nanmean(behaveMat);
y = nanstd(behaveMat);
zscores = (behaveMat - x)./y;
correlationMatrix = corr(zscores,'rows','pairwise');

imagesc(correlationMatrix)
set(gca,'YTickLabel',cellColNames(behaveStatsIdx),'YTick',1:length(behaveStatsIdx),...
    'XTickLabel',cellColNames(behaveStatsIdx),'XTick',1:length(behaveStatsIdx),...
    'XTickLabelRotation',45,'GridLineStyle','-','GridColor','black','GridAlpha',1,...
    'FontWeight','bold','LineWidth',1)

pbaspect([1 1 1])
skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0],1:256);
colormap(skalafell)
colorbar
caxis([-1 1]);

savefig('shiftedBehaviors')


%%


%%

% obtain transition 0'd kernels for each fly

% time of all transitions
accumulatedTime = [0; cumsum(lightDat(1:(end-1),1))];
transDiff = [0; diff(lightDat(:,2))];

% obtains light and dark assignments and time since transition stamps for every turn
for ii = 1:length(allDataStruct.tTime)
    turns = allDataStruct.tSequence{ii};
    time = allDataStruct.tTime{ii};
    %     wallDist = allDataStruct.tWallDist{ii};
    
    timeRel = zeros(length(time),1);
    for jj = 1:length(turns)
        % returns difference in time from current turn and last dark to light transition
        ldT = time(jj) - accumulatedTime(transDiff==1);
        % returns difference in time from current turn and last light to dark transition
        dlT = time(jj) - accumulatedTime(transDiff==-1);
        
        if ~isempty(min([min(dlT(dlT >= 0)) min(ldT(ldT >= 0))]))
            timeRel(jj) = min([min(dlT(dlT >= 0)) min(ldT(ldT >= 0))]);
            
        else
            timeRel(jj) = time(jj);
        end
    end
    
    allDataStruct.timeRel{ii} = timeRel;
    
    ldVect = logical(ldVect);
    timeRel = timeRel;
    if mod(ii,500) == 0
        [repmat('|',1,floor(ii/length(allDataStruct.tTime)*100)),...
            repmat(' ',1,100-floor(ii/length(allDataStruct.tTime)*100))]
    end
    
    
    % low temp light bins
    idx = allDataStruct.ldVect{ii} & lowTIdx{ii};
    turns = allDataStruct.tSequence{ii}(idx);
    normTime = timeRel(idx);
    time = allDataStruct.tTime{ii}(idx);
    walls = allDataStruct.tWallDist{ii}(idx);
    
    if length(turns)<=20
        continue
    end
    
    out = individualStatistics(allDataStruct);
    
    allDataStruct.lowLightTB{ii} = out.TB;
    allDataStruct.lowLightNumTurns{ii} = out.numTurns;
    allDataStruct.lowLightTBBSE{ii} = out.bse;
    allDataStruct.lowLightKernel{ii} = out.kernel;
    allDataStruct.lowLightKernelTurns{ii} = out.kernelTurns;
    allDataStruct.lowLightAuto{ii} = out.autoCorrelation;
    allDataStruct.lowLightAutoCI{ii} = out.autoCorrelationCI;
    allDataStruct.lowLightClump{ii} = out.clump;
    allDataStruct.lowLightSwitch(ii) = out.switch;
    allDataStruct.lowLightTBSplit(ii,:) = out.tbSplit;
    %     allDataStruct.lowLightWallPred(ii) = out.wallPred;
    %
    % high temp light bins
    idx = allDataStruct.ldVect{ii} & highTIdx{ii};
    turns = allDataStruct.tSequence{ii}(idx);
    normTime = allDataStruct.timeRel{ii}(idx);
    time = allDataStruct.timeRel{ii}(idx);
    %     walls = allDataStruct.tWallDist{ii}(idx);
    
    if length(turns)>10
        out = individualStatistics(turns,normTime,time);
    end
    
    %     assign data
    allDataStruct.highLightTB{ii} = out.TB;
    allDataStruct.highLightNumTurns{ii} = out.numTurns;
    allDataStruct.highLightTBBSE{ii} = out.bse;
    allDataStruct.highLightKernel{ii} = out.kernel;
    allDataStruct.highLightKernelTurns{ii} = out.kernelTurns;
    allDataStruct.highLightAuto{ii} = out.autoCorrelation;
    allDataStruct.highLightAutoCI{ii} = out.autoCorrelationCI;
    allDataStruct.highLightClump{ii} = out.clump;
    allDataStruct.highLightSwitch(ii) = out.switch;
    allDataStruct.highLightTBSplit(ii,:) = out.tbSplit;
    %     allDataStruct.highLightWallPred(ii) = out.wallPred;
    
    
    % low temp dark bins
    idx = ~allDataStruct.ldVect{ii} & lowTIdx{ii};
    turns = allDataStruct.tSequence{ii}(idx);
    normTime = allDataStruct.timeRel{ii}(idx);
    time = allDataStruct.timeRel{ii}(idx);
    %         walls = allDataStruct.tWallDist{ii}(idx);
    
    if length(turns)>10
        out = individualStatistics(turns,normTime,time);
    end
    
    %     assign data
    allDataStruct.lowDarkTB{ii} = out.TB;
    allDataStruct.lowDarkNumTurns{ii} = out.numTurns;
    allDataStruct.lowDarkTBBSE{ii} = out.bse;
    allDataStruct.lowDarkKernel{ii} = out.kernel;
    allDataStruct.lowDarkKernelTurns{ii} = out.kernelTurns;
    allDataStruct.lowDarkAuto{ii} = out.autoCorrelation;
    allDataStruct.lowDarkAutoCI{ii} = out.autoCorrelationCI;
    allDataStruct.lowDarkClump{ii} = out.clump;
    allDataStruct.lowDarkSwitch(ii) = out.switch;
    allDataStruct.lowDarkTBSplit(ii,:) = out.tbSplit;
    %     allDataStruct.lowDarkWallPred(ii) = out.wallPred;
    
    
    % high temp dark bins
    idx = ~allDataStruct.ldVect{ii} & highTIdx{ii};
    turns = allDataStruct.tSequence{ii}(idx);
    normTime = allDataStruct.timeRel{ii}(idx);
    time = allDataStruct.timeRel{ii}(idx);
    %         walls = allDataStruct.tWallDist{ii}(idx);
    
    if length(turns)>10
        out = individualStatistics(turns,normTime,time);
    end
    
    %     assign data
    allDataStruct.highDarkTB{ii} = out.TB;
    allDataStruct.highDarkNumTurns{ii} = out.numTurns;
    allDataStruct.highDarkTBBSE{ii} = out.bse;
    allDataStruct.highDarkKernel{ii} = out.kernel;
    allDataStruct.highDarkKernelTurns{ii} = out.kernelTurns;
    allDataStruct.highDarkAuto{ii} = out.autoCorrelation;
    allDataStruct.highDarkAutoCI{ii} = out.autoCorrelationCI;
    allDataStruct.highDarkClump{ii} = out.clump;
    allDataStruct.highDarkSwitch(ii) = out.switch;
    allDataStruct.highDarkTBSplit(ii,:) = out.tbSplit;
    %     allDataStruct.highDarkWallPred(ii) = out.wallPred;
    
end


%% Stats for all data
clf

% ldm compared to null
tb1 = cat(1,allDataStruct.highLightTB{:});
tb2 = cat(1,allDataStruct.highDarkTB{:});
ldm =  tb1 - tb2;
err1 = cat(1,allDataStruct.highLightTBBSE{:});
err2 = cat(1,allDataStruct.highLightTBBSE{:});

nullDist = normrnd(tb1,err1) - normrnd(tb1,err2);

figure
hold on
histogram(ldm,-0.5:0.005:0.5,'EdgeColor','none','Normalization','probability')
histogram(nullDist,-0.5:0.005:0.5,'EdgeColor','none','Normalization','probability')
xlim([-0.5,0.5])


figure
hold on
histogram(tb1,'EdgeColor','none','Normalization','probability')
histogram(tb2,'EdgeColor','none','Normalization','probability')
xlim([0,1])

shg

%% Generate Indexing variable

rubinAnno = importsplit('/Users/kyobikakaria/Desktop/Data/split_screen/split_Rubin_Annotation.csv');

genos = rubinAnno.Properties.VariableNames(2:end);
neurons = rubinAnno(:,1);
effectors = {'ISO','SHI','TRP'};

% this obtains the index of every occurance of a particular genotype separated by neuron type
genoIdx = cell(0);
for ii = 1:size(rubinAnno,1)
    logical1 = zeros(size(allDataStruct.tLabels,1),1);
    idx = logical(table2array(rubinAnno(ii,2:end)));
    g = genos(idx);
    for kk = 1:length(g)
        logical1(:,kk) = strcmp(g(kk),allDataStruct.tLabels(:,1));
    end
    genoIdx{ii} = logical1;
end

% this obtains the index of every occurance of a particular effector
logical2 = zeros(size(allDataStruct.tLabels,1),length(effectors));
for ii = 1:length(effectors)
    logical2(:,ii) = strcmp(effectors(ii),upper(allDataStruct.tLabels(:,2)));
end

logical2 = logical(logical2);



%% Create Network view of each geno/effector

reSamp = 100;

effectNames(1:3) = {'lightTBMADES','darkTBMADES','LDMES'};

effectNames(4:7) = {'lightTurnsES','darkTurnsES','lightTurnsMADES','darkTurnsMADES'};

effectNames(8:9) = {'lightTBSplitES','darkTBSplitES'};

effectNames(10:13) = {'lightTBSwitchES','darkTBSwitchES','lightTurnsMADES','darkTBSwitchMADES'};

effectNames(14:17) = {'lightTBClumpES','darkTBClumpES','lightTBClumpMADES','darkTBClumpMADES'};

for ii = 1:length(genoIdx)
    
    % index control or shibire
    controlIdx = any(genoIdx{ii},2) & logical2(:,1);
    shibireIdx = any(genoIdx{ii},2) & logical2(:,2);
    
    if sum(controlIdx) == 0 || sum(shibireIdx) == 0
        continue
    end
    
    % extract turn biases for each condition
    controlTBs = cat(1,allDataStruct.highLightTB{controlIdx});
    shibireTBs = cat(1,allDataStruct.highLightTB{shibireIdx});
    
    % extract dark turn biases for each condition
    controlTBsD = cat(1,allDataStruct.highDarkTB{controlIdx});
    shibireTBsD = cat(1,allDataStruct.highDarkTB{shibireIdx});
    
    lightTBMADES = (mad(shibireTBs,1) - mad(controlTBs,1))./...
        mad(controlTBs,1);
    
    darkTBMADES = (mad(shibireTBsD,1) - mad(controlTBsD,1))./...
        mad(controlTBsD,1);
    
    LDMES = (mad(shibireTBs - shibireTBsD,1) - mad(controlTBs - controlTBsD,1))./...
        mad(controlTBs - controlTBsD,1);
    
    effects(ii,1:3) = [lightTBMADES' darkTBMADES' LDMES'];
    
    
    
    %     % calculate mean-absolute deviation of turn bias distribution
    %     controlVar = nanstd(controlTBs);
    %     shibireVar = nanstd(shibireTBs);
    %
    %     % calculate SEM of MADs
    %     temp = nan(length(controlTBs),1);
    %     temp2 = nan(length(controlTBs),1);
    %     for jj = 1:reSamp
    %     temp(jj) = nanvar(randsample(controlTBs,length(controlTBs),1),1);
    %     temp2(jj) = nanvar(randsample(shibireTBs,length(shibireTBs),1),1);
    %     end
    %     controlVarSEM = nanstd(temp);
    %     shibireVarSEM = nanstd(temp2);
    %
    %     % calculate difference between control and shibire
    %     deltaVar(ii) = (shibireVar - controlVar)./controlVar;
    %     deltaVarSEM(ii) = sqrt(shibireVarSEM.^2 + controlVarSEM.^2)./controlVar;
    %
    %     % calculate ldm of turn biases in the high temperature
    %     % I am now defining ldm as the variance of the differences distribution
    %
    %     % calculate ldm
    %     controlDiffs = controlTBs - controlTBsD;
    %     controlLDM = nanvar(controlDiffs);
    %     nullVarC = nan(100,1);
    %     for jj = 1:100
    %         nullDist = normrnd(controlTBs,controlTBsErrs) -  normrnd(controlTBs,controlTBsDErrs);
    %         nullVarC(jj) = nanvar(nullDist);
    %     end
    %     controlResVar = controlLDM-mean(nullVarC);
    %
    %     shibireDiffs = shibireTBs - shibireTBsD;
    %     shibireLDM = nanvar(shibireDiffs);
    %     nullVarS = nan(100,1);
    %     for jj = 1:100
    %         nullDist = normrnd(shibireTBs,shibireTBsErrs) -  normrnd(shibireTBs,shibireTBsDErrs);
    %         nullVarS(jj) = nanvar(nullDist);
    %     end
    %     shibireResVar = shibireLDM - mean(nullVarS);
    %
    %
    %     ldmShiftRe = nan(reSamp,1);
    %     for jj = 1:reSamp
    %         conReSamp = randsample(controlDiffs,length(controlDiffs),1);
    %         controlResVarReSamp = nanvar(conReSamp) - mean(nullVarC);
    %         shiReSamp = randsample(shibireDiffs,length(controlDiffs),1);
    %         shibireResVarReSamp = nanvar(shiReSamp) - mean(nullVarS);
    %         ldmShiftRe(jj) = log2(shibireResVarReSamp./controlResVarReSamp);
    %     end
    %
    %     % calculate fold change from control to shibire
    %     ldmShift(ii) = log2(shibireResVar./controlResVar);
    %     ldmShiftE(ii) = nanstd(ldmShiftRe);
    %
    % calculate activity variables
    
    % extract turn numbers for each condition
    controlTurns = cat(1,allDataStruct.highLightNumTurns{controlIdx});
    shibireTurns = cat(1,allDataStruct.highLightNumTurns{shibireIdx});
    
    % extract dark turn numbers for each condition
    controlTurnsD = cat(1,allDataStruct.highDarkNumTurns{controlIdx});
    shibireTurnsD = cat(1,allDataStruct.highDarkNumTurns{shibireIdx});
    
    
    lightTurnsES = (nanmean(shibireTurns) - nanmean(controlTurns))./...
        nanmean(controlTurns);
    
    darkTurnsES = (nanmean(shibireTurnsD) - nanmean(controlTurns))./...
        nanmean(controlTurns);
    
    lightTurnsMADES = (mad(shibireTurnsD,1) - mad(controlTurnsD,1))./...
        mad(controlTurnsD,1);
    
    darkTurnsMADES = (mad(shibireTurnsD,1) - mad(controlTurnsD,1))./...
        mad(controlTurnsD,1);
    
    effects(ii,4:7) = [lightTurnsES',darkTurnsES',lightTurnsMADES',darkTurnsMADES'];
    
    
    
    
    %     % calculate difference between control and shibire
    %     deltaTurnVar(ii) = (shibireVar - controlVar)./controlVar;
    %     deltaTurnSEM(ii) = sqrt(shibireVarSEM.^2 + controlVarSEM.^2)./controlVar;
    %
    
    
    % intra-individual variation calculation
    
    % extract turn numbers for each condition
    controlTBSplit = cat(1,allDataStruct.highLightTBSplit(controlIdx,:));
    shibireTBSplit = cat(1,allDataStruct.highLightTBSplit(shibireIdx,:));
    
    % extract turn numbers for each condition
    controlTBSplitD = cat(1,allDataStruct.highDarkTBSplit(controlIdx,:));
    shibireTBSplitD = cat(1,allDataStruct.highDarkTBSplit(shibireIdx,:));
    
    lightTBSplitES = (nanstd(diff(shibireTBSplit,[],2)) - nanstd(diff(controlTBSplit,[],2)))./...
        nanstd(diff(controlTBSplit,[],2));
    
    darkTBSplitES = (nanstd(diff(shibireTBSplitD,[],2)) - nanstd(diff(controlTBSplitD,[],2)))./...
        nanstd(diff(controlTBSplitD,[],2));
    
    effects(ii,8:9) = [lightTBSplitES',darkTBSplitES'];
    
    
    
    % switchiness calculation
    
    % extract turn numbers for each condition
    controlSwitch = cat(1,allDataStruct.highLightSwitch(controlIdx));
    shibireSwitch = cat(1,allDataStruct.highLightSwitch(shibireIdx));
    
    % extract turn numbers for each condition
    controlSwitchD = cat(1,allDataStruct.highDarkSwitch(controlIdx));
    shibireSwitchD = cat(1,allDataStruct.highDarkSwitch(shibireIdx));
    
    lightTBSwitchES = (nanmean(shibireSwitch) - nanmean(controlSwitch))./...
        nanstd(controlSwitch);
    
    darkTBSwitchES = (nanmean(shibireSwitchD) - nanmean(controlSwitchD))./...
        nanstd(controlSwitchD);
    
    lightTBSwitchMADES = (mad(shibireSwitch,1) - mad(controlSwitch,1))./...
        mad(controlSwitch,1);
    
    darkTBSwitchMADES = (mad(shibireSwitchD,1) - mad(controlSwitchD,1))./...
        mad(controlSwitchD,1);
    
    effects(ii,10:13) = [lightTBSwitchES',darkTBSwitchES',lightTBSwitchMADES',darkTBSwitchMADES'];
    
    
    
    % clumpiness calculation
    
    % extract turn numbers for each condition
    controlClump = cat(1,allDataStruct.highLightClump{controlIdx});
    shibireClump = cat(1,allDataStruct.highLightClump{shibireIdx});
    
    % extract turn numbers for each condition
    controlClumpD = cat(1,allDataStruct.highDarkClump{controlIdx});
    shibireClumpD = cat(1,allDataStruct.highDarkClump{shibireIdx});
    
    lightTBClumpES = (nanmean(shibireClump) - nanmean(controlClump))./...
        nanstd(controlSwitch);
    
    darkTBClumpES = (nanmean(shibireClumpD) - nanmean(controlClumpD))./...
        nanstd(controlClumpD);
    
    lightTBClumpMADES = (mad(shibireClump,1) - mad(controlClump,1))./...
        mad(controlSwitch,1);
    
    darkTBClumpMADES = (mad(shibireClumpD,1) - mad(controlClumpD,1))./...
        mad(controlClumpD,1);
    
    effects(ii,14:17) = [lightTBClumpES',darkTBClumpES',lightTBClumpMADES',darkTBClumpMADES'];
    
    
    
end

%% call network drawing function

subplot(2,3,4)
drawingNetwork(LDMES)
title('LDM')
subplot(2,3,5)
drawingNetwork(deltaVar,deltaVarSEM)
title('TB Variance')
subplot(2,3,6)
drawingNetwork(deltaTurnVar,deltaTurnSEM)
title('Delta Activity')

%%

preds = [allDataStruct.lowDarkWallPred' allDataStruct.lowLightWallPred'];



preds = [allDataStruct.highDarkWallPred' allDataStruct.highLightWallPred'];
