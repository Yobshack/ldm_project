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
% for ii = 1:length(files)
% turnDirs{ii} = extractTurns(files{ii},lightDat)
% end

%% Make large cell array to handle all data that I can query for groups
% 
% lightTimes = turnDirs{1}.tTime{1}(turnDirs{1}.tSequence{1});
% darkTimes = turnDirs{1}.tTime{1}(~turnDirs{1}.tSequence{1});
% 
% clf
% hold on
% colormap(parula)
% h = histogram(diff(lightTimes),0:0.2:20,'FaceColor','blue','EdgeColor','blue',...
%     'Normalization','probability');
% h.EdgeAlpha = 1;
% h = histogram(diff(darkTimes),0:0.2:20,'FaceColor','red','EdgeColor','red',...
%     'Normalization','probability');
% h.EdgeAlpha = 0.2;
% 
% plot(0.1:0.2:20,poisspdf(1:100,poissfit(h.BinCounts(1:100))))
shg
%% Make large cell array to handle all data that I can query for groups

files = dir('*_processed.mat');
tmpC = cell(0);
for ii = 1:length(files)
    load(files(ii).name)
    
    % need to heal a few of the labels files because the number of roi's was not 120 on some trays
    if size(turnDirs.tLabels,1) > size(turnDirs.tSequence,1)
        turnDirs.tLabels = turnDirs.tLabels(1:size(turnDirs.tSequence,1),:);
    elseif size(turnDirs.tLabels,1) < size(turnDirs.tSequence,1)
        turnDirs.tLabels(end:size(turnDirs.tSequence,1),:) = turnDirs.tLabels(end,:);
    end     
    
    dataCell = cell(size(turnDirs.tSequence,1),12);
    for kk = 1:size(turnDirs.tSequence,1)
        
        % separate geno and effector
        labs = turnDirs.tLabels(kk,:);
        labs(6) = cell(1);
        labsSplit = strsplit(char(labs(1)),'_');
        labs = [labsSplit labs(2:(end-1))];
        dataCell(kk,:) = [turnDirs.tSequence(kk),turnDirs.tTime(kk),labs,...
            {turnDirs.dateAndTime},turnDirs.avgWallDist(kk),...
            logical(turnDirs.avgLight(kk)),turnDirs.proArea(kk)];
    end
    tmpC{ii} = dataCell;
end

cellColNames = {'turnSequence','turnTimes','genotype','effector','sex',...
    'treatment','box','maze','dateAndTime','wallDistance','lightStatus',...
    'totalMazeArea'};

% concatenate all data into the same cell array
allDataCell = cat(1,tmpC{:});


% filter for proper maze size
a = cell2mat(allDataCell(:,12));
idx = a > (median(a) - 2*mad(a)) & a < (median(a) + 2*mad(a));

allDataCell = allDataCell(idx,:);

% time in minutes
lowT = [0 60]*60;
highT = [120 180]*60;

% turnBias calculation
for ii = 1:size(allDataCell,1)
    
    % return index for all low temp turns
    idxL = allDataCell{ii,2} >= lowT(1) & allDataCell{ii,2} < lowT(2);
    % total turns low
    allDataCell(ii,13) = num2cell(length(allDataCell{ii,1}(idxL)));
    % light turns low
    allDataCell(ii,14) = num2cell(length(allDataCell{ii,1}(idxL & allDataCell{ii,11})));
    % dark turns low
    allDataCell(ii,15) = num2cell(length(allDataCell{ii,1}(idxL & ~allDataCell{ii,11})));

    idxH = allDataCell{ii,2} >= highT(1) & allDataCell{ii,2} < highT(2);
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
    turns = allDataCell{ii,1}(idxL);
    iti = diff(allDataCell{ii,2}(idxL));
    allDataCell(ii,25) = num2cell(mad(iti,1)/(time/length(turns)));
    % light clump low
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
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,31} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    % light switch low
    turns = allDataCell{ii,1}(idxL & allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,32} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    % dark switch low
    turns = allDataCell{ii,1}(idxL & ~allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,33} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    
    % total switch low
    turns = allDataCell{ii,1}(idxH);
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,34} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    % light switch low
    turns = allDataCell{ii,1}(idxH & allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,35} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    % dark switch low
    turns = allDataCell{ii,1}(idxH & ~allDataCell{ii,11});
    div = 2*length(turns)*nanmean(turns)*(1-nanmean(turns));
    allDataCell{ii,36} = sum((turns(1:end-1)+turns(2:end))==1)/div;
    
    % wall following idx
    
    % total wall low
    walls = allDataCell{ii,10}(idxL);
    allDataCell{ii,37} = nanmean(walls);
    % light wall low
    walls = allDataCell{ii,10}(idxL & allDataCell{ii,11});
    allDataCell{ii,38} = nanmean(walls);
    % dark wall low
    walls = allDataCell{ii,10}(idxL & ~allDataCell{ii,11});
    allDataCell{ii,39} = nanmean(walls);
    
    % wall following idx
    
    turns = allDataCell{ii,1};
    % total wall high left
    walls = allDataCell{ii,10}(idxH & ~turns);
    allDataCell{ii,40} = nanmean(walls);
    % total wall high left
    walls = allDataCell{ii,10}(idxH & turns);
    allDataCell{ii,41} = nanmean(walls);
    % light wall high left
    walls = allDataCell{ii,10}(idxH & ~turns & allDataCell{ii,11});
    allDataCell{ii,42} = nanmean(walls);
    % light wall high right
    walls = allDataCell{ii,10}(idxH & turns & allDataCell{ii,11});
    allDataCell{ii,43} = nanmean(walls);
    % dark wall high left
    walls = allDataCell{ii,10}(idxH & ~turns & ~allDataCell{ii,11});
    allDataCell{ii,44} = nanmean(walls);
    % dark wall high right
    walls = allDataCell{ii,10}(idxH & turns & ~allDataCell{ii,11});
    allDataCell{ii,45} = nanmean(walls);
    
end
% get indexes for low and high temperature times
for ii = 1:length(allDataStruct.tTime)
    lowTIdx{ii} = allDataStruct.tTime{ii} > lowT(1)*60 & allDataStruct.tTime{ii} <= lowT(2)*60;
    highTIdx{ii} = allDataStruct.tTime{ii} > highT(1)*60 & allDataStruct.tTime{ii} <= highT(2)*60;
end

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
