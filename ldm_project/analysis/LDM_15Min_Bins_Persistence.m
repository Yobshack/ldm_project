%% Kyobi Skutt-Kakaria
% 05.02.2019 - Created


%% Select files to analyze

% opens a ui to select mat files to analyze
fullPath = pwd;
files = uigetfile('*.mat','MultiSelect','on');
if iscell(files) ~= 1
    files = {files};
end
files = sort(files);

% loads an experimental design file, script expects a n x 3 array where n is the number of bins in
% the stimulus sequence. col 1 is the bin length in seconds, col 2 is a logical vector indicating
% whether the `lights were on or off in that bin. col 3 is a vector of pwm values that reflect the
% pwm setting of the arduino controlling the led light board.
lightDat = dlmread('lightSequenceScreen_15Mins.txt');
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

files = uigetfile('*.mat','MultiSelect','on');
files = sort(files);
[allDataCell, cellColNames] = parseProcessedFiles(files);

%% Calculate behavioral metrics


for ii = 1:size(allDataCell,1)
    idx = allDataCell{ii,11} & allDataCell{ii,2} < 120*60;
    allDataCell{ii,13} = turnbias(allDataCell{ii,1}(idx));
    
    idx = allDataCell{ii,11} & allDataCell{ii,2} > 120*60;
    allDataCell{ii,14} = turnbias(allDataCell{ii,1}(idx));
    
    idx = ~allDataCell{ii,11} & allDataCell{ii,2} < 120*60;
    allDataCell{ii,15} = turnbias(allDataCell{ii,1}(idx));
    
    idx = ~allDataCell{ii,11} & allDataCell{ii,2} > 120*60;
    allDataCell{ii,16} = turnbias(allDataCell{ii,1}(idx));
end

date = cellfun(@(x) cell2mat(x(1:2)), allDataCell(:,9),'Un',0);
[unique_dates,~,date] = unique(date);
[unique_maze,~,maze] = unique(cell2mat(allDataCell(:,8)));
tbLight1 = cell2mat(allDataCell(:,13));
tbLight2 = cell2mat(allDataCell(:,14));
tbDark1 = cell2mat(allDataCell(:,15));
tbDark2 = cell2mat(allDataCell(:,16));

outMat = nan(1,4,6);
c = 1;
for ii = 1:length(unique_maze)
    for jj = 1:length(unique_dates)
        idx = maze == ii & date == jj;
        try
        outMat(ii,jj,:) = cat(3,tbLight1(idx,1),tbLight2(idx,1),...
            tbDark1(idx,1),tbDark2(idx,1),...
            tbLight1(idx,1) - tbDark1(idx,1),...
            tbLight2(idx,1) - tbDark2(idx,1)); 
        end
        c = c+1;
    end

end

outMat(outMat==0) = nan;

clear lightCurve darkCurve ldmCurve
numRe = 10000;
for ii = 1:numRe
    idx = randi(size(outMat,1),size(outMat,1),1);
    tmp = corr(outMat(idx,:,1),outMat(idx,:,2),'type','Spearman','rows','pairwise');
    lightCurve(ii,:) = tmp(1,:);
    
    tmp = corr(outMat(idx,:,3),outMat(idx,:,4),'type','Spearman','rows','pairwise');
    darkCurve(ii,:) = tmp(1,:);
    
    tmp = corr(outMat(idx,:,5),outMat(idx,:,6),...
        'type','Spearman','rows','pairwise');
    ldmCurve(ii,:) = tmp(1,:);
end



col = [0 0.4 1; 0.8 0.2 0; 0 0 0];
hErr3 = errorbar(nanmedian(ldmCurve),nanstd(ldmCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(3,:),...
    'MarkerFaceColor',col(3,:));
hold on

hErr2 = errorbar(nanmedian(darkCurve),nanstd(darkCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(2,:),...
    'MarkerFaceColor',col(2,:));

hErr1 = errorbar(nanmedian(lightCurve),nanstd(lightCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(1,:),...
    'MarkerFaceColor',col(1,:));

hold off

xlim([0 5])
ylim([0 1])
xticks([1:4])
set([hErr1,hErr2,hErr3],'MarkerSize',20)
legend('LDM','Dark','Light')

shg

n = size(outMat,1);

day1LDM = outMat(:,1,6);
day2LDM = outMat(:,2,6);

eles = {'ldm_correlation_resamples','dark_correlation_resamples','light_correlation_resamples',...
    'median_correlations','day1_ldm','day2_ldm','n'};
LDM_Figure1i = {ldmCurve,darkCurve,lightCurve,[nanmedian(ldmCurve),...
    nanmedian(darkCurve),nanmedian(lightCurve)],day1LDM,day2LDM,n,eles};

save('LDM_Figure1i.mat','LDM_Figure1i');

%%

for ii = 1:size(allDataCell,1)
    idx = allDataCell{ii,11} & allDataCell{ii,2} < 120*60;
    allDataCell{ii,13} = wallDist(allDataCell{ii,10}(idx),allDataCell{ii,1}(idx));
    
    idx = allDataCell{ii,11} & allDataCell{ii,2} > 120*60;
    allDataCell{ii,14} = wallDist(allDataCell{ii,10}(idx),allDataCell{ii,1}(idx));
    
    idx = ~allDataCell{ii,11} & allDataCell{ii,2} < 120*60;
    allDataCell{ii,15} = wallDist(allDataCell{ii,10}(idx),allDataCell{ii,1}(idx));
    
    idx = ~allDataCell{ii,11} & allDataCell{ii,2} > 120*60;
    allDataCell{ii,16} = wallDist(allDataCell{ii,10}(idx),allDataCell{ii,1}(idx));
end


date = cellfun(@(x) cell2mat(x(1:2)), allDataCell(:,9),'Un',0);
[unique_dates,~,date] = unique(date);
[unique_maze,~,maze] = unique(cell2mat(allDataCell(:,8)));
tbLight1 = cell2mat(allDataCell(:,13));
tbLight2 = cell2mat(allDataCell(:,14));
tbDark1 = cell2mat(allDataCell(:,15));
tbDark2 = cell2mat(allDataCell(:,16));

outMat = nan(1,4,6);
c = 1;
for ii = 1:length(unique_maze)
    for jj = 1:length(unique_dates)
        idx = maze == ii & date == jj;
        try
        outMat(ii,jj,:) = cat(3,tbLight1(idx,1),tbLight2(idx,1),...
            tbDark1(idx,1),tbDark2(idx,1),...
            tbLight1(idx,1) - tbDark1(idx,1),...
            tbLight2(idx,1) - tbDark2(idx,1)); 
        end
        c = c+1;
    end

end

outMat(outMat==0) = nan;

clear lightCurve darkCurve ldmCurve
numRe = 10000;
for ii = 1:numRe
    idx = randi(size(outMat,1),size(outMat,1),1);
    tmp = corr(outMat(idx,:,1),outMat(idx,:,2),'type','Spearman','rows','pairwise');
    lightCurve(ii,:) = tmp(1,:);
    
    tmp = corr(outMat(idx,:,3),outMat(idx,:,4),'type','Spearman','rows','pairwise');
    darkCurve(ii,:) = tmp(1,:);
    
    tmp = corr(outMat(idx,:,5),outMat(idx,:,6),...
        'type','Spearman','rows','pairwise');
    ldmCurve(ii,:) = tmp(1,:);
end



col = [0 0.4 1; 0.8 0.2 0; 0 0 0];
hErr3 = errorbar(nanmedian(ldmCurve),nanstd(ldmCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(3,:),...
    'MarkerFaceColor',col(3,:));
hold on

hErr2 = errorbar(nanmedian(darkCurve),nanstd(darkCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(2,:),...
    'MarkerFaceColor',col(2,:));

hErr1 = errorbar(nanmedian(lightCurve),nanstd(lightCurve),'CapSize',0,'LineWidth',2,'Marker','o','Color',col(1,:),...
    'MarkerFaceColor',col(1,:));

hold off

xlim([0 5])
ylim([0 1])
xticks([1:4])
set([hErr1,hErr2,hErr3],'MarkerSize',20)
legend('LDM','Dark','Light')

shg

n = size(outMat,1);

day1LDM = outMat(:,1,6);
day2LDM = outMat(:,2,6);

scatter(day1LDM,day2LDM);

eles = {'ldm_correlation_resamples','dark_correlation_resamples','light_correlation_resamples',...
    'median_correlations','day1_ldm','day2_ldm','n'};
LDM_FigureS1j = {lightCurve,darkCurve,ldmCurve,[nanmedian(ldmCurve),...
    nanmedian(darkCurve),nanmedian(lightCurve)],day1LDM,day2LDM,n,eles};

save('LDM_FigureS1j','LDM_FigureS1j');



  