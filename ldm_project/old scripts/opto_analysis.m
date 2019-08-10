%% Kyobi Skutt-Kakaria
% 02.27.2018 - Created
% 02.27.2018 - updated


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
lightDat = dlmread('lightSequence_Opto90Secs.txt');
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

files = dir('*_processed.mat');
[allDataCell, cellColNames] = parseProcessedFiles(files);

%% Calculate behavioral metrics
% Takes a while to run this on a large data set
% order is light - low, dark - low, light - high, dark - high
% if ~isempty(dir('check1.mat'))
%     load('check1.mat')
% else

% factors - time bins, temperature bins, light condition, animals
tempTimes = {[0 360*60]};
    %cellColNames(13:14) = {'LightLow','DarkLow','LightHigh','DarkHigh'};
    parfor hh = 1:size(allDataCell,1)
        hh
        data = allDataCell(hh,:);
        bins = unique(lightDat(:,3));
        [dataMat{hh},errorMat{hh}] = generateBehavioralMetrics3(data,lightDat,tempTimes,bins);
    end

    array = cell(1,14);
    for ii = 1:size(allDataCell,1)          
        data = dataMat{ii};
        error = errorMat{ii};
        array(ii,:) = [allDataCell(ii,:) data error];
    end
% end

array(:,3) = upper(array(:,3));

data = cat(1,array{:,13});
error = cat(1,array{:,14});



%% Generate Rank VBE plot

idx1 = strcmp(array(:,6),'-ATR') & strcmp(array(:,3),'SS02252');
idx2 = strcmp(array(:,6),'+ATR') & strcmp(array(:,3),'SS02252');
bin = 1;
%%
generateDistPlots([data(idx1,:,1,:,bin),data(idx1,:,2,:,bin)],[error(idx1,:,1,:,bin) error(idx1,:,2,:,bin)]);
generateDistPlots([data(idx2,:,1,:,bin),data(idx2,:,2,:,bin)],[error(idx2,:,1,:,bin) error(idx2,:,2,:,bin)]);
%%

generateDiffPlot2([data(idx1,:,1,:,bin),data(idx1,:,2,:,bin)],[error(idx1,:,1,:,bin) error(idx1,:,2,:,bin)])
generateDiffPlot2([data(idx2,:,1,:,bin),data(idx2,:,2,:,bin)],[error(idx2,:,1,:,bin) error(idx2,:,2,:,bin)])

%% Generate MI plot
light = 0;
figure
hold on
colors = {[0 0 0],[0.8 0.2 0],[0 0.6 0.2]};
bins = unique(lightDat(:,3));
bin = 1;
    idx = idx1; %ugi == ii & x;% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        hTurn1 = generateTurnTMIPlot2(array(idx,:),lightDat,light,colors{1},bins(bin));
        val = (max(hTurn1.mainLine.YData)-min(hTurn1.mainLine.YData))/2+min(hTurn1.mainLine.YData);
        [~,ix] = min(abs(hTurn1.mainLine.YData(hTurn1.mainLine.XData <= 0) - val));
        line([hTurn1.mainLine.XData(1), hTurn1.mainLine.XData(ix), hTurn1.mainLine.XData(ix),...
            hTurn1.mainLine.XData(ix)],...
            [hTurn1.mainLine.YData(ix), hTurn1.mainLine.YData(ix), hTurn1.mainLine.YData(ix), 0],...
            'Color',colors{1})
        drawnow
    end
    idx = idx2; %ugi == ii & x% & cell2mat((array(:,7))) == kk;
    if sum(idx) > 0
        hTurn2 = generateTurnTMIPlot2(array(idx,:),lightDat,light,colors{2},bins(bin));
        val = (max(hTurn2.mainLine.YData)-min(hTurn2.mainLine.YData))/2+min(hTurn2.mainLine.YData);
        [~,ix] = min(abs(hTurn2.mainLine.YData(hTurn2.mainLine.XData <= 0) - val));
        line([hTurn2.mainLine.XData(1), hTurn2.mainLine.XData(ix), hTurn2.mainLine.XData(ix),...
            hTurn2.mainLine.XData(ix)],...
            [hTurn2.mainLine.YData(ix), hTurn2.mainLine.YData(ix), hTurn2.mainLine.YData(ix), 0],...
            'Color',colors{2})
        drawnow
    end
    
%%

%figure
colors = {[0 0 0],[0.8 0.2 0],[0 0.6 0.2]}; 
alpha = 0.5;

figure
hold on

behave = 2;
light = 1;

bin = 3

f = @(p,x) p(1) + ((p(2)-p(1)) ./ ((1 + p(3)*exp(-p(4)*x)).^(1/p(5))));
h = fittype(@(p1,p2,p3,p4,p5,x) p1 + ((p2-p1) ./ ((1 + p3*exp(-p4*x)).^(1/p5))));

idx = idx1;
subArray = array(idx,:);
hError1 = generateTTAPlot2(subArray,lightDat,behave,colors{1},light,bin);
% ind = hError1.XData > 0;
% xVect = hError1.XData(1):0.1:hError1.XData(end);
% a1 = fit(hError1.XData',hError1.YData',h);
% hPlot1 = plot(xVect,a1(xVect),'Color',colors{1},'LineWidth',2);
% hLine1 = line([a1(3), a1(3), -5, a1(3)],[0, f(a1,a1(3)), f(a1,a1(3)), f(a1,a1(3))],...
%     'Color',[colors{1} alpha],'LineWidth',2,'LineStyle','-','LineJoin','round','LineWidth',0.5);
% 

drawnow

idx = idx2
subArray = array(idx,:);
hError2 = generateTTAPlot2(subArray,lightDat,behave,colors{2},light,bin);
% a2 = fit(hError2.XData',hError2.YData',h);
% xVect = hError2.XData(1):0.1:hError2.XData(end);
% hPlot2 = plot(xVect,a2(xVect),'Color',colors{2},'LineWidth',2);
% hLine2 = line([a2(3), a2(3), -5, a2(3)],[0, f(a2,a2(3)), f(a2,a2(3)), f(a2,a2(3))],...
%     'Color',[colors{2} alpha],'LineWidth',2,'LineStyle','-','LineJoin','round','LineWidth',0.5);


drawnow