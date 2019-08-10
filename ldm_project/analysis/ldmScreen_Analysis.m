%% LDM Manuscript Screen Analysis
% 12.06.2018
% Kyobi Skutt-Kakaria

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

%% Neuron annotation file

% import rubin annotation
filepath = '/Volumes/LaCie/Data/split_screen/split_Gal4/split_Rubin_Annotation.csv';
importedData = importdata(filepath);

names = importedData.textdata(2:end,1);
genonames = importedData.textdata(1,2:end);


shorterNames = {'P-FN_M/P','P-FN_A','P-FN_D','P-FN_V','P-EN','P-EG','E-PG','E-PG_T','P-F-G_S','P-F-R',...
    'PF-LC','LPs-P',...
    'P_6-8-P_9','Delta7','Sps-P','IbSps-P'};

names(1:length(shorterNames)) = shorterNames;

% image of neurons by drivers included in screen
sub = importedData.data(:,[1:47 81:91]);
id = nan(size(sub,2),1);
for ii = 1:size(sub,2)
    try
        a = find(sub(:,ii));
        id(ii) = a(1);
    end
end
[~,sortI] = sort(id);
hImage = imagesc(sub(:,sortI));

colormap(gray)
shg


%% Make large cell array to handle all data that I can query for groups

files = dir('*_processed.mat');
[allDataCell, cellColNames] = parseProcessedFiles(files);


%% Calculate behavioral metrics
% Takes a while to run this on a large data set
% order is light - low, dark - low, light - high, dark - high

times = {[0 60]*60,[120 180]*60};
turnFilter = 50;
bins = 1;
allData = cell(0); allError = allData;
for hh = 1:size(allDataCell,1)
    hh
    [allData{hh}, allError{hh}] = generateBehavioralMetrics3(allDataCell(hh,:),lightDat,times,...
        bins,turnFilter);
    hh
end

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

array(:,3) = upper(array(:,3));

save('LDM_Analysis.mat','array')

clear

%% Load saved analysis .mat file

load('LDM_Analysis.mat')


%% Figure 1n_Supps

% Sets neuron index
allLines = cell(0);
for kk = 1:23
    try
        genoN = kk;
        
        % Sets behavior index 1 = activity, 2 = turn bias, 3 = switchiness,
        % 4 = wall following, 5 = clumpiness
        behave = 2;
        
        % retrieve effectors
        effectors = array(:,4);
        
        % retrieve all genotypes that are part of the neuron label
        ngn = genonames(logical(importedData.data(genoN,:)));
        if kk <= 16
            ngn = ngn(~contains(ngn,'R'));
        end
        
        % time block 1 = low, 2 = high
        time = 2;
        
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10;
        
        ldmCorrected = nan(length(ngn),2);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        for ii = 1:length(ngn)
            
            try
                
                % Set index for specific experiment for shibire and control
                idx1 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff1);
                idx2 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff2);
                
                flyBehavior1 = cat(1,array{idx1,13});
                flyBehavior2 = cat(1,array{idx2,13});
                
                % get individual experiment id's
                [experiments1,~,experimentInd1] = uniqueRowsCA(cat(1,array{idx1,9}));
                [experiments2,~,experimentInd2] = uniqueRowsCA(cat(1,array{idx2,9}));
                
                % Set turn bias to variable and zscore turn biases
                normVals1 = [nanmean(flyBehavior1(:,2,1,time)),nanstd(flyBehavior1(:,2,1,time))];
                turnBiasInLight = (flyBehavior1(:,2,1,time) - normVals1(1))./normVals1(2);
                normVals2 = [nanmean(flyBehavior1(:,2,2,time)),nanstd(flyBehavior1(:,2,2,time))];
                turnBiasInDark = (flyBehavior1(:,2,2,time) - normVals2(1))./normVals2(2);
                
                % Null Model = All animals are p = 0.5
                turnsInLight = round(flyBehavior1(:,1,1,time)*(times{2}(2)-times{2}(1))/2);  % turns in 30 mins;
                turnsInDark = round(flyBehavior1(:,1,2,time)*(times{2}(2)-times{2}(1))/2);
                nullLDM = nan(numRe,1);
                for jj = 1:numRe
                    nullLight = binornd(turnsInLight,0.5)./turnsInLight;
                    nullLight = (nullLight - 0.5)./normVals1(2);
                    nullDark = binornd(turnsInDark,0.5)./turnsInDark;
                    nullDark = (nullDark - 0.5)./normVals2(2);
                    nullLDM(jj) = nanvar(nullLight-nullDark);
                end
                
                % Calculate LDM by taking variance of LDM and subtracting null variance
                ldmDistribution = turnBiasInLight - turnBiasInDark;
                ldmCorrected(ii,1) = sqrt(nanvar(ldmDistribution) - nanmean(nullLDM));
                ldmCorrectedE{ii,1} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nanmean(nullLDM)),...
                    ldmDistribution);
                
                % Set turn bias to variable
                normVals1 = [nanmean(flyBehavior2(:,2,1,time)),nanstd(flyBehavior2(:,2,1,time))];
                turnBiasInLight = (flyBehavior2(:,2,1,time) - normVals1(1))./normVals1(2);
                normVals2 = [nanmean(flyBehavior2(:,2,2,time)),nanstd(flyBehavior2(:,2,2,time))];
                turnBiasInDark = (flyBehavior2(:,2,2,time) - normVals2(1))./normVals2(2);
                
                % Null Model = All animals are p = 0.5
                turnsInLight = round(flyBehavior2(:,1,1,time)*(times{2}(2)-times{2}(1))/2);  % turns in 30 mins;
                turnsInDark = round(flyBehavior2(:,1,2,time)*(times{2}(2)-times{2}(1))/2);
                nullLDM = nan(numRe,1);
                for jj = 1:numRe
                    nullLight = binornd(turnsInLight,0.5)./turnsInLight;
                    nullLight = (nullLight - 0.5)./normVals1(2);
                    nullDark = binornd(turnsInDark,0.5)./turnsInDark;
                    nullDark = (nullDark - 0.5)./normVals2(2);
                    nullLDM(jj) = nanvar(nullLight-nullDark);
                end
                
                % Calculate LDM by taking variance of LDM and subtracting null variance
                ldmDistribution = turnBiasInLight - turnBiasInDark;
                ldmCorrected(ii,2) = sqrt(nanvar(ldmDistribution) - nanmean(nullLDM));
                ldmCorrectedE{ii,2} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nanmean(nullLDM)),...
                    ldmDistribution);
                
                ldmDelta = (ldmCorrected(:,1) - ldmCorrected(:,2))./ldmCorrected(:,2);
                for jj = 1:size(ldmCorrectedE,1)
                    ldmDeltaDist{jj} = (ldmCorrectedE{jj,1} - ldmCorrectedE{jj,2})...
                        ./ldmCorrectedE{jj,2};
                end
                
                ii
            end
            
        end
        
        
        % Make figures
        
        close all
        data = cat(2,ldmDeltaDist{:});
        idx = ~all(isnan(data),1);
        data = data(:,idx);
        
        figure
        line([0 size(data,2)+1],[0 0],'color',[0.5 0.5 0.5])
        hold on
        hViolin = violin(data,[],'facecolor',[0.3 0.3 0.3],'facealpha',0.5);
        
        ylabel('LDM')
        title(names{kk})
        
        pVals = sum(data>0)./size(data,1);
        pValsString = cell(6,1);
        pValsString(1:length(pVals)) = {''};
        
        pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
        pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
        pValsString(pVals < 0.001) = {'***'};
        
        for ll = 1:length(pVals)
            text([ll ll],[0.5 0.5],pValsString{ll},'HorizontalAlignment','center','FontSize',36)
        end
        
        
        set(gca,'XTick',1:size(data,2),'XTickLabel',ngn(idx),'FontSize',24,...
            'YLim',[-1 1],'XLim',[0 size(data,2)+1],'XTickLabelRotation',45)
        
        pbaspect([0.5 + (size(data,2))/5,1,1])
        hFigure = gcf;
        hFigure.Position =  [680   558   560   420];
        drawnow
        print(strcat('LDM_figureS1n_',num2str(kk),'.pdf'),'-dpdf');
        
        
        allLines{kk} = data;
        
    end
end


save('LDM_figureS1n.mat','allLines','names')

%% Make draft of Figure 1N
% Sets neuron index
allLines = cell(0);
for kk = 1:23
    try
        genoN = kk;
        
        % Sets behavior index 1 = activity, 2 = turn bias, 3 = switchiness,
        % 4 = wall following, 5 = clumpiness
        behave = 2;
        
        % retrieve effectors
        effectors = array(:,4);
        
        % retrieve all genotypes that are part of the neuron label
        ngn = genonames(logical(importedData.data(genoN,:)));
        if kk <= 16
            ngn = ngn(~contains(ngn,'R'));
        end
        
        % time block 1 = low, 2 = high
        time = 2;
        
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10000;
        
        ldmCorrected = nan(length(ngn),2);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        
        try
            % Set index for specific experiment for shibire and control
            idx1 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff1);
            idx2 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff2);
            
            flyBehavior1 = cat(1,array{idx1,13});
            flyBehavior2 = cat(1,array{idx2,13});
            flyError1 = cat(1,array{idx1,14});
            flyError2 = cat(1,array{idx2,14});
            
            % get individual experiment id's
            [experiments1,~,experimentInd1] = uniqueRowsCA(cat(1,array{idx1,9}));
            [experiments2,~,experimentInd2] = uniqueRowsCA(cat(1,array{idx2,9}));
            
            % Set turn bias to variable and zscore turn biases
            normVals1 = [nanmean(flyBehavior1(:,2,1,time)),nanstd(flyBehavior1(:,2,1,time))];
            turnBiasInLight = (flyBehavior1(:,2,1,time) - normVals1(1))./normVals1(2);
            turnBiasInLightError = (flyError1(:,2,2,time))./normVals1(2);
            normVals2 = [nanmean(flyBehavior1(:,2,2,time)),nanstd(flyBehavior1(:,2,2,time))];
            turnBiasInDark = (flyBehavior1(:,2,2,time) - normVals2(1))./normVals2(2);
            turnBiasInDarkError = (flyError1(:,2,2,time))./normVals2(2);
            
            %                 % Null Model = All animals are p = 0.5
            %                 turnsInLight = round(flyBehavior1(:,1,1,time)*(times{2}(2)-times{2}(1))/2);  % turns in 30 mins;
            %                 turnsInDark = round(flyBehavior1(:,1,2,time)*(times{2}(2)-times{2}(1))/2);
            %                 nullLDM = nan(numRe,1);
            %                 for jj = 1:numRe
            %                     nullLight = binornd(turnsInLight,0.5)./turnsInLight;
            %                     nullLight = (nullLight - 0.5)./normVals1(2);
            %                     nullDark = binornd(turnsInDark,0.5)./turnsInDark;
            %                     nullDark = (nullDark - 0.5)./normVals2(2);
            %                     nullLDM(jj) = nanvar(nullLight-nullDark);
            %                 end
            %
            % Subtract expected standard error
            nullLDMErr = nanmean(turnBiasInLightError.^2 + turnBiasInDarkError.^2);
            
            % Calculate LDM by taking variance of LDM and subtracting null variance
            ldmDistribution = turnBiasInLight - turnBiasInDark;
            ldmCorrected(kk,1) = sqrt(nanvar(ldmDistribution) - nullLDMErr);
            ldmCorrectedE{kk,1} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nullLDMErr),...
                ldmDistribution);
            
            %                 % Calculate LDM by taking variance of LDM and subtracting null variance
            %                 ldmDistribution = turnBiasInLight - turnBiasInDark;
            %                 ldmCorrected(kk,1) = sqrt(nanvar(ldmDistribution) - nanmean(nullLDM));
            %                 ldmCorrectedE{kk,1} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nanmean(nullLDM)),...
            %                     ldmDistribution);
            %
            
            % Set turn bias to variable
            normVals1 = [nanmean(flyBehavior2(:,2,1,time)),nanstd(flyBehavior2(:,2,1,time))];
            turnBiasInLight = (flyBehavior2(:,2,1,time) - normVals1(1))./normVals1(2);
            turnBiasInLightError = (flyError2(:,2,2,time))./normVals1(2);
            normVals2 = [nanmean(flyBehavior2(:,2,2,time)),nanstd(flyBehavior2(:,2,2,time))];
            turnBiasInDark = (flyBehavior2(:,2,2,time) - normVals2(1))./normVals2(2);
            turnBiasInDarkError = (flyError2(:,2,2,time))./normVals2(2);
            
            %                 % Null Model = All animals are p = 0.5
            %                 turnsInLight = round(flyBehavior2(:,1,1,time)*(times{2}(2)-times{2}(1))/2);  % turns in 30 mins;
            %                 turnsInDark = round(flyBehavior2(:,1,2,time)*(times{2}(2)-times{2}(1))/2);
            %                 nullLDM = nan(numRe,1);
            %                 for jj = 1:numRe
            %                     nullLight = binornd(turnsInLight,0.5)./turnsInLight;
            %                     nullLight = (nullLight - 0.5)./normVals1(2);
            %                     nullDark = binornd(turnsInDark,0.5)./turnsInDark;
            %                     nullDark = (nullDark - 0.5)./normVals2(2);
            %                     nullLDM(jj) = nanvar(nullLight-nullDark);
            %                 end
            %
            %                 % Calculate LDM by taking variance of LDM and subtracting null variance
            %                 ldmDistribution = turnBiasInLight - turnBiasInDark;
            %                 ldmCorrected(kk,2) = sqrt(nanvar(ldmDistribution) - nanmean(nullLDM));
            %                 ldmCorrectedE{kk,2} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nanmean(nullLDM)),...
            %                     ldmDistribution);
            %
            % Subtract expected standard error
            nullLDMErr = nanmean(turnBiasInLightError.^2 + turnBiasInDarkError.^2);
            
            % Calculate LDM by taking variance of LDM and subtracting null variance
            ldmDistribution = turnBiasInLight - turnBiasInDark;
            ldmCorrected(kk,2) = sqrt(nanvar(ldmDistribution) - nullLDMErr);
            ldmCorrectedE{kk,2} = bootstrp(numRe,@(x) sqrt(nanvar(x) - nullLDMErr),...
                ldmDistribution);
            
            
            ldmDelta = (ldmCorrected(:,1) - ldmCorrected(:,2))./ldmCorrected(:,2);
            for jj = 1:size(ldmCorrectedE,1)
                ldmDeltaDist{jj} = (ldmCorrectedE{jj,1} - ldmCorrectedE{jj,2})...
                    ./ldmCorrectedE{jj,2};
            end
            
            kk
        end
        
        
        data = cat(2,ldmDeltaDist{:});
        idx = ~all(isnan(data));
        data = data(:,idx);
        
        allLines{kk} = data;
        
    end
end

data = cell(0);
for ii = 1:length(allLines)
    data{ii} = allLines{ii}(:);
    pVals(ii) = sum(data{ii}>0)./size(data{ii},1);
end


figure
line([0 size(data,2)+1],[0 0],'color',[0.5 0.5 0.5])
hold on

violin(data)
set(gca,'XTickLabel',names)


ylabel('LDM')



pValsString = cell(length(data),1);
pValsString(1:length(pVals)) = {''};

pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
pValsString(pVals < 0.001) = {'***'};

for ll = 1:length(pVals)
    text([ll ll],[0.5 0.5],pValsString{ll},'HorizontalAlignment','center','FontSize',16)
end


set(gca,'XTick',1:size(data,2),'FontSize',12,...
    'YLim',[-1 1],'XLim',[0 size(data,2)+1],'XTickLabelRotation',45)

pbaspect([0.5 + (size(data,2))/5,1,1])
hFigure = gcf;
hFigure.Position =  [680   558   1600   600];
drawnow
data{:,24} = cellfun(@(x) nanmedian(x),data(:,1:23));
data{:,25} = names;
print('-bestfit',strcat('LDM_figure1n','.pdf'),'-dpdf');
save('LDM_figure1n.mat','data')


%% Figure S2e All lines


% Sets neuron index
allLines = cell(0);
scores = cell(0);
for kk = 1:23
    %%
    try
        genoN = kk;
        
        % Sets behavior index 1 = activity, 2 = turn bias, 3 = switchiness,
        % 4 = wall following, 5 = clumpiness
        behave = 2;
        
        % retrieve effectors
        effectors = array(:,4);
        
        % retrieve all genotypes that are part of the neuron label
        ngn = genonames(logical(importedData.data(genoN,:)));
        if kk <= 16
            ngn = ngn(~contains(ngn,'R'));
        end
        
        neuronNames{kk} = ngn;
        
        % time block 1 = low, 2 = high
        time = 2;
        
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10;
        
        ldmCorrected = nan(length(ngn),2);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        
        for ii = 1:length(ngn)
            
            try
                % Set index for specific experiment for shibire and control
                idx1 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff1);
                idx2 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff2);
                
                % Subset by these indexes
                flyBehavior1 = cat(1,array{idx1,13});
                activityI = flyBehavior1(:,1,2,1) > 0.03;
                flyBehavior1 = flyBehavior1(activityI,:,:,:);
                flyBehavior2 = cat(1,array{idx2,13});
                activityI = flyBehavior2(:,1,2,1) > 0.03;
                flyBehavior2 = flyBehavior2(activityI,:,:,:);
                
                % Set turn bias to variable and zscore turn biases for high temperature
                time = 2;
                normVals1 = [nanmean(flyBehavior1(:,2,1,time)),nanstd(flyBehavior1(:,2,1,time))];
                turnBiasInLightHigh = (flyBehavior1(:,2,1,time) - normVals1(1))./normVals1(2);
                %                 turnBiasInLightHigh = flyBehavior1(:,2,1,time);
                normVals2 = [nanmean(flyBehavior1(:,2,2,time)),nanstd(flyBehavior1(:,2,2,time))];
                turnBiasInDarkHigh = (flyBehavior1(:,2,2,time) - normVals2(1))./normVals2(2);
                %                 turnBiasInDarkHigh = flyBehavior1(:,2,2,time);
                
                % Set turn bias to variable and zscore turn biases for low temperature
                time = 1;
                normVals1 = [nanmean(flyBehavior1(:,2,1,time)),nanstd(flyBehavior1(:,2,1,time))];
                turnBiasInLightLow = (flyBehavior1(:,2,1,time) - normVals1(1))./normVals1(2);
                %                 turnBiasInLightLow = flyBehavior1(:,2,1,time);
                normVals2 = [nanmean(flyBehavior1(:,2,2,time)),nanstd(flyBehavior1(:,2,2,time))];
                turnBiasInDarkLow = (flyBehavior1(:,2,2,time) - normVals2(1))./normVals2(2);
                %                 turnBiasInDarkLow = flyBehavior1(:,2,2,time);
                % %
                % %                 % Similarity will be calculated by calculating the correlation coefficient between
                % %                 % the zscored values and their pre-temperature matched controls
                
                flyBiases1 = [turnBiasInLightLow, turnBiasInDarkLow,...
                    turnBiasInLightHigh, turnBiasInDarkHigh];
                
                %                 [~,ind] = sort(abs(turnBiasInLightLow - turnBiasInDarkLow));
                %                 ind = ind(length(ind)/2:length(ind));
                %                 flyBiases1 = flyBiases1(ind,:);
                
                % Control
                
                % Set turn bias to variable and zscore turn biases for high temperature
                time = 2;
                normVals1 = [nanmean(flyBehavior2(:,2,1,time)),nanstd(flyBehavior2(:,2,1,time))];
                turnBiasInLightHigh = (flyBehavior2(:,2,1,time) - normVals1(1))./normVals1(2);
                %                 turnBiasInLightHigh = flyBehavior2(:,2,1,time);
                normVals2 = [nanmean(flyBehavior2(:,2,2,time)),nanstd(flyBehavior2(:,2,2,time))];
                turnBiasInDarkHigh = (flyBehavior2(:,2,2,time) - normVals2(1))./normVals2(2);
                %                 turnBiasInDarkHigh = flyBehavior2(:,2,2,time);
                
                % Set turn bias to variable and zscore turn biases for low temperature
                time = 1;
                normVals1 = [nanmean(flyBehavior2(:,2,1,time)),nanstd(flyBehavior2(:,2,1,time))];
                turnBiasInLightLow = (flyBehavior2(:,2,1,time) - normVals1(1))./normVals1(2);
                %                 turnBiasInLightLow = flyBehavior2(:,2,1,time);
                normVals2 = [nanmean(flyBehavior2(:,2,2,time)),nanstd(flyBehavior2(:,2,2,time))];
                turnBiasInDarkLow = (flyBehavior2(:,2,2,time) - normVals2(1))./normVals2(2);
                %                 turnBiasInDarkLow = flyBehavior2(:,2,2,time);
                
                flyBiases2 = [turnBiasInLightLow, turnBiasInDarkLow,...
                    turnBiasInLightHigh, turnBiasInDarkHigh];
                
                %                 [~,ind] = sort(abs(turnBiasInLightLow - turnBiasInDarkLow));
                %                 ind = ind(length(ind)/2:length(ind));
                %                 flyBiases2 = flyBiases2(ind,:);
                
                scores{kk}{ii} = nan(numRe,4);
                for hh = 1:numRe
                    indexes = randi(size(flyBiases1,1),size(flyBiases1,1),1);
                    corrmap = corr(flyBiases1(indexes,:),'rows','pairwise','type','Spearman');
                    corrs1 = [corrmap(3,1), corrmap(3,2)]; 
                    
                    indexes = randi(size(flyBiases2,1),size(flyBiases2,1),1);
                    corrmap = corr(flyBiases2(indexes,:),'rows','pairwise','type','Spearman');
                    corrs2 = [corrmap(3,1) - corrmap(3,2), corrmap(4,1) - corrmap(4,2)];
                    
                    scores{kk}{ii}(hh,:) = [corrs1,corrs2];
                    
                    %                     indexes = randi(size(flyBiases1,1),size(flyBiases1,1),1);
                    %                     [minimums,minI] = min(abs(flyBiases1(indexes,3) - flyBiases1(indexes,1:2)),[],2);
                    %                     scores{ii}(hh,1) = (nansum(minI==1)./numel(minI)' - 0.5)*2;
                    %                     [minimums,minI] = min(abs(flyBiases1(indexes,4) - flyBiases1(indexes,1:2)),[],2);
                    %                     scores{ii}(hh,2) = (nansum(minI==1)./numel(minI)' - 0.5)*2;
                    %                     indexes = randi(size(flyBiases2,1),size(flyBiases2,1),1);
                    %                     [minimums,minI] = min(abs(flyBiases2(indexes,3) - flyBiases2(indexes,1:2)),[],2);
                    %                     scores{ii}(hh,3) = (nansum(minI==1)./numel(minI)' - 0.5)*2;
                    %                     [minimums,minI] = min(abs(flyBiases2(indexes,4) - flyBiases2(indexes,1:2)),[],2);
                    %                     scores{ii}(hh,4) = (nansum(minI==1)./numel(minI)' - 0.5)*2;
                end
                
                
            end
            
        end
        
    end
end
%%
% Make figures

close all
total = 0;
clf
for kk = 1:length(scores)
    kk
    dataLight = cell(0);
    for ii = 1:length(scores{kk})
        dataLight{ii} = (scores{kk}{ii}(:,1) - scores{kk}{ii}(:,3));
    end
    
    dataLight = cat(2,dataLight{:});
    
    pVals = sum(dataLight>0)./size(dataLight,1);
    pValsString = cell(size(dataLight,2),1);
    pValsString(1:length(pVals)) = {''};
    
    pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
    pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
    pValsString(pVals < 0.001) = {'***'};
    
    for ii = 1:size(dataLight,2)
        subplot(6,10,ii+total)
        
        line([0 2],[0 0],'color',[0.5 0.5 0.5])
        hold on
        hViolin = violin(dataLight(:,ii),[],'facecolor',[0 0.8 1],'facealpha',1);
        set(gca,'XTick',1,'XTickLabel',neuronNames{kk}(ii),'FontSize',8,...
            'YLim',[-1 1],'XLim',[0.5 1.5],'XTickLabelRotation',45)
        
        pbaspect([1 1 1])
        hFigure = gcf;
        ylabel('LDM Direction')
        title(names{kk})
        
        text([1 1],[0.25 0.25],pValsString{ii},'HorizontalAlignment',...
            'center','FontSize',6)
        
        total = total+1;
        drawnow
    end
end
%%
print(strcat('LDM_figureS2e_',num2str(kk),'.pdf'),'-dpdf','-bestfit');

close all

dataDark = cell(0);
for ii = 1:length(scores)
    dataDark{ii} = (scores{ii}(:,2) - scores{ii}(:,4));
end

dataDark = cat(2,dataDark{:});

figure
line([0 size(dataDark,2)+1],[0 0],'color',[0.5 0.5 0.5])
hold on

hViolin = violin(dataDark,[],'facecolor',[1 0.1 0],'facealpha',1,'medc',...
    'b');

set(gca,'XTick',1:size(dataDark,2),'XTickLabel',ngn,'FontSize',24,...
    'YLim',[-1 1],'XLim',[0 size(dataDark,2)+1],'XTickLabelRotation',45)

pbaspect([0.5 + (size(dataDark,2))/5,1,1])
hFigure = gcf;
ylabel('LDM Direction')
title(names{kk})

pVals = sum(dataDark<0)./size(dataDark,1);
pValsString = cell(size(dataDark,2),1);
pValsString(1:length(pVals)) = {''};

pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
pValsString(pVals < 0.001) = {'***'};

for ll = 1:length(pVals)
    text([ll ll],[0.75 0.75],pValsString{ll},'HorizontalAlignment','center','FontSize',36)
end


pbaspect([0.5 + (size(dataDark,2))/5,1,1])
hFigure = gcf;
hFigure.Position =  [680   558   560   420];
drawnow
print(strcat('LDM_figureS2e_',num2str(kk+23),'.pdf'),'-dpdf','-bestfit');


allLines(kk,:) = {dataLight,dataDark};




%% Figure 2e Summary


% Sets neuron index
allLines = cell(0);
allBehave = cell(0);
         
actFilter = 0.03; % turns per second
flyN = nan(23,2);

for kk = 1:23
    %%
    
    try
        genoN = kk;
        
        % Sets behavior index 1 = activity, 2 = turn bias, 3 = switchiness,
        % 4 = wall following, 5 = clumpiness
        behave = 2;
        
        % retrieve effectors
        effectors = array(:,4);
        
        % retrieve all genotypes that are part of the neuron label
        ngn = genonames(logical(importedData.data(genoN,:)));
        if kk <= 16
            ngn = ngn(~contains(ngn,'R'));
        end
        
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10000;
        
        ldmCorrected = nan(length(ngn),2);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        scores = cell(0);
        %         for ii = 1:length(ngn)
        ii = 1;
        
        try
            % Set index for specific experiment for shibire and control
            idx1 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff1);
            idx2 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff2);
            
            % Zscore turn biases for each genotype
            flyBehavior1 = cat(1,array{idx1,13});
            flyError1 = cat(1,array{idx1,14});
            activityI = flyBehavior1(:,1,2,1) > actFilter;
            flyBehavior1 = flyBehavior1(activityI,:,:,:);
            
            flyBehavior2 = cat(1,array{idx2,13});
            flyError2 = cat(1,array{idx2,14});
            activityI = flyBehavior2(:,1,2,1) > actFilter;
            flyBehavior2 = flyBehavior2(activityI,:,:,:);
            
            
            % Experiment
            
            % Set biases to variable
            flyBiases1 = reshape(flyBehavior1(:,2,:,:),size(flyBehavior1,1),4); % Low Light, Low Dark, High Light, High Dark
            flyBiasesE1 = reshape(flyError1(:,2,:,:),size(flyError1,1),4);
            
            noNans = all(~isnan(flyBiases1),2);
            flyBiases1 = flyBiases1(noNans,:);
            flyBiasesE1 = flyBiasesE1(noNans,:);
            
            
            % Control
            
            % Set biases to variable
            flyBiases2 = reshape(flyBehavior2(:,2,:,:),size(flyBehavior2,1),4); % Low Light, Low Dark, High Light, High Dark
            flyBiasesE2 = reshape(flyError2(:,2,:,:),size(flyError2,1),4); % Low Light, Low Dark, High Light, High Dark
            
            noNans = all(~isnan(flyBiases2),2);
            flyBiases2 = flyBiases2(noNans,:);
            flyBiasesE2 = flyBiasesE2(noNans,:);
            flyN(kk,1:2) = [size(flyBiases1,1),size(flyBiases2,1)] ;
            
            scores{ii} = nan(numRe,8);
            for hh = 1:numRe
                indexes = randi(size(flyBiases1,1),size(flyBiases1,1),1);
                corrmap = corr(flyBiases1(indexes,:),'type','Spearman');
                fCorr = @(x,y) [x,y];
                corrs1 = [fCorr(corrmap(3,1),corrmap(3,2)), fCorr(corrmap(4,2),corrmap(4,1))];
                
                indexes = randi(size(flyBiases2,1),size(flyBiases2,1),1);
                corrmap = corr(flyBiases2(indexes,:),'type','Spearman');
                corrs2 = [fCorr(corrmap(3,1),corrmap(3,2)), fCorr(corrmap(4,2),corrmap(4,1))];
                
                scores{ii}(hh,:) = [corrs1,corrs2];
                
            end
            kk
        end
        
        
        % Extract data
        
        f = @(x,y) [x(:,1:4)./nanmean(nanmean((x(:,[1,3])),2)),...
            y(:,1:4)./nanmean(nanmean((y(:,[1,3])),2))] ;
        outData = f(scores{ii}(:,1:4),scores{ii}(:,5:8));
        
        dataLight = cell(0);
        for ii = 1:length(scores)
            dataLight{ii} = outData(:,2) - outData(:,6);
        end        
        dataLight = cat(2,dataLight{:});      
           
        dataDark = cell(0);
        for ii = 1:length(scores)
            dataDark{ii} = outData(:,4) - outData(:,8);
        end        
        dataDark = cat(2,dataDark{:});
               
        allLines(kk,:) = {dataLight,dataDark};
        allBehave{kk} = {flyBiases1,flyBiases2};
        
        if kk == 20
            r1 = scores;
        elseif kk == 7
            epg = scores;
        elseif kk == 1
            pfn = scores;
        elseif kk == 11
            plc = scores;
            
        end
        
    end
end

r11 = r1{1}(:,[1 5 2 6 3 7 4 8]);
epg1 = epg{1}(:,[1 5 2 6 3 7 4 8]);
pfn1 = pfn{1}(:,[1 5 2 6 3 7 4 8]);
plc1 = plc{1}(:,[1 5 2 6 3 7 4 8]);

labs = {'LL - Shi','LL - ctl', ' LD - shi', 'LD - ctl','DD - shi','DD - ctl','DL - shi','DL - ctl'};
eles = {'Data as shown','geno','n','labels'};

violin(r11,'BoxStyle','filled')
xticklabels(labs)
ylabel('Correlation Coefficient')
shg
mat = {r11,names(20),flyN(20,:),eles};
save('LDM_Figure_3X_rh1','mat')

violin(epg1,'BoxStyle','filled')
xticklabels(labs)
ylabel('Correlation Coefficient')
mat = {epg1,names(7),flyN(7,:),eles};
save('LDM_Figure_3X_epg','mat')

violin(pfn1,'BoxStyle','filled')
xticklabels(labs)
ylabel('Correlation Coefficient')
mat = {pfn1,names(1),flyN(1,:),eles};
save('LDM_Figure_3X_pfn','mat')

violin(plc1,'BoxStyle','filled')
xticklabels(labs)
ylabel('Correlation Coefficient')
mat = {plc1,names(11),flyN(11,:),eles};
save('LDM_Figure_3X_pflc','mat')


%%
dataLight = cell(0);

f = @(x) [sum(x>0),sum(x<0)]./size(x,1);
for ii = 1:length(allLines)
    dataLight{ii} = allLines{ii,1}(:);
    sums = f(dataLight{ii});
    pVals(ii) = min(sums);
end

idx = ~cellfun(@(x) isempty(x),dataLight);

dataLight = dataLight(idx);

figure
line([0 size(dataLight,2)+1],[0 0],'color',[0.5 0.5 0.5])
hold on

violin(dataLight,[],'facecolor',[0 0.8 1],'facealpha',1)
set(gca,'XTickLabel',names(idx))


ylabel('Role in Light Bias')

pValsString = cell(length(dataLight),1);
pValsString(1:length(pVals)) = {''};

pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
pValsString(pVals < 0.001) = {'***'};

for ll = 1:length(pVals)
    text([ll ll],[0.5 0.5],pValsString{ll},'HorizontalAlignment','center','FontSize',16)
end


set(gca,'XTick',1:size(dataLight,2),'FontSize',16,...
    'YLim',[-0.5 1],'XLim',[0 size(dataLight,2)+1],'XTickLabelRotation',45)

pbaspect([0.5 + (size(dataLight,2))/5,1,1])
hFigure = gcf;
hFigure.Position =  [680   558   1600   600];
drawnow
filename = strcat('LDM_figureS2e_47');
print('-fillpage',strcat(filename,'.pdf'),'-dpdf');
save(strcat(filename,'.mat'),'dataLight','names')


f = @(x) [sum(x>0),sum(x<0)]./size(x,1);
dataDark = cell(0);
for ii = 1:length(allLines)
    dataDark{ii} = allLines{ii,2}(:);
    sums = f(dataDark{ii});
    pVals(ii) = min(sums);
end


dataDark = dataDark(idx);

figure
line([0 size(dataDark,2)+1],[0 0],'color',[0.5 0.5 0.5])
hold on

violin(dataDark,[],'facecolor',[1 0.1 0],'facealpha',1,'medc',...
    'b')
set(gca,'XTickLabel',names(idx))


ylabel('Role in Dark Bias')

pValsString = cell(length(dataDark),1);
pValsString(1:length(pVals)) = {''};

pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
pValsString(pVals < 0.001) = {'***'};

for ll = 1:length(pVals)
    text([ll ll],[0.5 0.5],pValsString{ll},'HorizontalAlignment','center','FontSize',16)
end


set(gca,'XTick',1:size(dataDark,2),'FontSize',16,...
    'YLim',[-0.5 1],'XLim',[0 size(dataDark,2)+1],'XTickLabelRotation',45)

pbaspect([0.5 + (size(dataDark,2))/5,1,1])
hFigure = gcf;
hFigure.Position =  [680   558   1600   600];
drawnow
filename = strcat('LDM_figureS2e_48');
print('-fillpage',strcat(filename,'.pdf'),'-dpdf');
save(strcat(filename,'.mat'),'dataDark','names')

figure
x = cellfun(@(x) mean(x),dataDark);
y = cellfun(@(x) mean(x),dataLight);
scatter(x,y)
text(x,y,names)
line([-0.5 1],[0 0],'color',[0 0 0 0.25],'LineWidth',0.5)
line([0 0],[-0.5 1],'color',[0 0 0 0.25],'LineWidth',0.5)
set(gca,'YLim',[-0.25 0.75],'XLim',[-0.25,0.75])
pbaspect([1,1,1])
ylabel('Role In Light Bias')
xlabel('Role In Dark Bias')


shg

eles = {'Light resamples','Dark Resamples','Light Medians','Dark Medians',"N's",'genotypes'};
outMat = {dataLight,dataDark,cellfun(@nanmedian,dataLight),cellfun(@nanmedian,dataDark),flyN,names,eles};
save('LDM_Figure3b-e','outMat')




%% Make draft of Figure 2N
% Sets neuron index
allLines = cell(0);
f = @(x) (x - nanmean(x))./nanstd(x);
for kk = 1:23
    
    try
        genoN = kk;
        
        % Sets behavior index 1 = activity, 2 = turn bias, 3 = switchiness,
        % 4 = wall following, 5 = clumpiness
        behave = 2;
        
        % retrieve effectors
        effectors = array(:,4);
        
        % retrieve all genotypes that are part of the neuron label
        ngn = genonames(logical(importedData.data(genoN,:)));
        if kk <= 16
            ngn = ngn(~contains(ngn,'R'));
        end
        
        % light block 1 = light, 2 = dark
        light = 1;
        
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10;
        
        ldmCorrected = nan(length(ngn),4);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        
        try
            % Set index for specific experiment for shibire and control
            idx1 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff1);
            idx2 = ismember(upper(array(:,3)),ngn) & ismember(effectors,eff2);
            
            flyBehavior1 = cat(1,array{idx1,13}); % dims: 1 - fly, 2 - behavior, 3 - light (1 light 2 darkREFREF), 4 - temperature (low,high)
            flyBehavior2 = cat(1,array{idx2,13});
            
        end
        % Subset by these indexes
        flyBehavior1 = cat(1,array{idx1,13});
        activityI = flyBehavior1(:,1,2,1) > 0.025;
        flyBehavior1 = flyBehavior1(activityI,:,:,:);
        flyBehavior2 = cat(1,array{idx2,13});
        activityI = flyBehavior2(:,1,2,1) > 0.025;
        flyBehavior2 = flyBehavior2(activityI,:,:,:);
        
        
        % Experimental Light
        x = flyBehavior1(:,2,1,1);
        y = flyBehavior1(:,2,1,2);
        xx = flyBehavior1(:,1,1,1);
        yy = flyBehavior1(:,1,1,2);
        tdmCorrected(:,1) = bootstrp(numRe,@(x,y,xx,yy) ldm(x,y,xx,yy,1800,100),x,y,xx,yy);
        ldm(x,y,xx,yy,1800,100)
        
        % Experimental Dark
        x = flyBehavior1(:,2,2,1);
        y = flyBehavior1(:,2,2,2);
        xx = flyBehavior1(:,1,2,1);
        yy = flyBehavior1(:,1,2,2);
        tdmCorrected(:,2) = bootstrp(numRe,@(x,y,xx,yy) ldm(x,y,xx,yy,1800,100),x,y,xx,yy);
        
        % Control Light
        x = flyBehavior2(:,2,1,1);
        y = flyBehavior2(:,2,1,2);
        xx = flyBehavior2(:,1,1,1);
        yy = flyBehavior2(:,1,1,2);
        tdmCorrected(:,3) = bootstrp(numRe,@(x,y,xx,yy) ldm(x,y,xx,yy,1800,100),x,y,xx,yy);
        
        % Control Dark
        x = flyBehavior2(:,2,2,1);
        y = flyBehavior2(:,2,2,2);
        xx = flyBehavior2(:,1,2,1);
        yy = flyBehavior2(:,1,2,2);
        tdmCorrected(:,4) = bootstrp(numRe,@(x,y,xx,yy) ldm(x,y,xx,yy,1800,100),x,y,xx,yy);
        
        allLines{kk} = tdmCorrected;
        
    end
    
    figure
    violin(tdmCorrected)
    
    kk
    pause(3)
end
%%

data = cell(0);
for ii = 1:length(allLines)
    data{ii} = allLines{ii}(:);
    pVals(ii) = sum(data{ii}>0)./size(data{ii},1);
end


figure
line([0 size(data,2)+1],[0 0],'color',[0.5 0.5 0.5])
hold on

violin(data)
set(gca,'XTickLabel',names)


ylabel('LDM')



pValsString = cell(length(data),1);
pValsString(1:length(pVals)) = {''};

pValsString(pVals < 0.05 & pVals >= 0.01) = {'*'};
pValsString(pVals < 0.01 & pVals >= 0.001) = {'**'};
pValsString(pVals < 0.001) = {'***'};

for ll = 1:length(pVals)
    text([ll ll],[0.5 0.5],pValsString{ll},'HorizontalAlignment','center','FontSize',16)
end


set(gca,'XTick',1:size(data,2),'FontSize',12,...
    'YLim',[-1 1],'XLim',[0 size(data,2)+1],'XTickLabelRotation',45)

pbaspect([0.5 + (size(data,2))/5,1,1])
hFigure = gcf;
hFigure.Position =  [680   558   1600   600];
drawnow
data{:,24} = cellfun(@(x) nanmedian(x),data(:,1:23));
data{:,25} = names;
print('-bestfit',strcat('LDM_figure1n','.pdf'),'-dpdf');
save('LDM_figure1n.mat','data')


