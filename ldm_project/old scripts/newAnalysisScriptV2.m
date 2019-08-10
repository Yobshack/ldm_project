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

%% Make large cell array to handle all data that I can query for groups

files = dir('*_processed.mat');



%% Calculate behavioral metrics

% create choosable data parameters
temps = {'low','high'};
tempTimes = {[0 60]*60,[120 180]*60};
light = {'light','dark','total'};
turnFilter = 30;

% Generate cell array with behavioral metric for each individual to be indexed later
allDataCellArray = cell(1000,length(temps),length(light));
for hh = 1300:size(allDataCell,1)
    data = allDataCell(hh,:);



for ii = 1:length(temps)
    for jj = 1:length(light)      
        idx1 = data{2} > tempTimes{ii}(1) & data{2} <= tempTimes{ii}(2);
        switch jj
            case 1 
                idx2 = data{11};
            case 2
                idx2 = ~data{11};
            case 3
                idx2 = ones(length(idx1),1);
        end
        idx = idx1 & idx2;
        turns = data{1}(idx);
        walls = data{10}(idx);
        tStamps = data{2}(idx);
        if length(turns) > turnFilter

        allDataCellArray{hh,ii,jj} = {numturns(turns),turnbias(turns),...
            switchiness(turns), wallDist(walls,turns),clumpiness(tStamps,lightDat)};
        end
    end
end
hh
end


%% Generate simple figures
% this section is a nested loop to process each group for all group statistics
% NEED TO MAKE INTO A PROPER FUNCTION


tempCon = find(strcmp(temps,'high'));

datMat = nan(2,5,2);
for kk = 1:2
tempArray = cat(1,allDataCellArray{:,tempCon,kk});

for jj = 1:size(tempArray,2)
    m = cat(1,tempArray{:,jj});
    if size(m,2) == 1
        m1 = bootstrp(100,@nanmean,m);
        m2 = bootstrp(100,@(x) nanstd(x(:,1))/nanmean(x(:,1)),m);
    else
        % this treats every data point as a normal distribution with error set by earlier methods
        m1 = bootstrp(100,@(x) nanmean(normrnd(x(:,1),x(:,2))),m);
        f = @(x) nanstd(x)/nanmean(x);
        m2 = bootstrp(100,@(x) f(normrnd(x(:,1),x(:,2))),m);
    end
    datMat(kk,jj,1) = nanmean(m(:,1));
    datMat(kk,jj,2) = std(m1);
    datMat(kk,jj,3) = nanstd(m(:,1))/nanmean(m(:,1));
    datMat(kk,jj,4) = nanstd(m2);
    
    hist(m(:,1),30)
    drawnow
end
end

barwitherr(datMat(1:2,:,4)',datMat(1:2,:,3)')
shg
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
