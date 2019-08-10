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


%% Make draft of Figure 1N
% Sets neuron index
clear n groups

separate_by_driver = true;
eff1 = 'SHI'; % TRP for dTrpA1, % SHI for shibire
eff2 = 'ISO';
numRe = 10000;
ldmDeltaDist = cell(0);
ldm_experimental_groups = cell(0);
time = 2;
geno_names = cell(0);
bias_turns_experiment = cell(0);
bias_turns_control = bias_turns_experiment;

c = 1;

for kk = 1:23
    
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
        
        
        
              
        if separate_by_driver
            num_drivers = length(ngn);
        else
            num_drivers = 1;
        end
        
        for ii = 1:num_drivers
       
            
            
            % Set index for specific experiment for shibire and control
            if separate_by_driver
                idx1 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff1);
                idx2 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff2);
                 if sum(idx1) <= 50  || sum(idx2) <= 50
                     continue
                 end
            else
                idx1 = cell(0);
                idx2 = cell(0);
                for hh = 1:length(ngn)                   
                    idx1{hh} = ismember(upper(array(:,3)),ngn{hh}) & ismember(effectors,eff1);
                    idx2{hh} = ismember(upper(array(:,3)),ngn{hh}) & ismember(effectors,eff2);
                    
                    if sum(idx1{hh}) <= 50  || sum(idx2{hh}) <= 50 
                        idx1{hh} = false(length(idx1{hh}),1);
                        idx2{hh} = idx1{hh};
                    end                  
                end
                
                idx1 = any(cell2mat(idx1),2);
                idx2 = any(cell2mat(idx2),2);            
            end
            
            flyBehavior1 = cat(1,array{idx1,13});
            flyBehavior2 = cat(1,array{idx2,13});
            flyError1 = cat(1,array{idx1,14});
            flyError2 = cat(1,array{idx2,14});
            
            % Experimental basic stats
            bias = [flyBehavior1(:,2,1,time),...
                flyBehavior1(:,2,2,time)];
            activity = round([flyBehavior1(:,1,1,time) * 60 * 30,...
                flyBehavior1(:,1,2,time) * 60 * 30]); % needs to be manually edited if time bins are changed at the moment
            bias_turns_experiment{c,1} = [bias,activity];
            
       
            % Control basic stats
            bias = [flyBehavior2(:,2,1,time),...
                flyBehavior2(:,2,2,time)];
            activity = round([flyBehavior2(:,1,1,time) * 60 * 30,...
                flyBehavior2(:,1,2,time) * 60 * 30]); % needs to be manually edited if time bins are changed at the moment
            bias_turns_control{c,1} = [bias,activity];      
            
            
            % EXPERIMENTAL GROUP
            
            % Set turn bias to variable and zscore turn biases
            normVals1 = [nanmean(flyBehavior1(:,2,1,time)),nanstd(flyBehavior1(:,2,1,time))];
            turnBiasInLight = (flyBehavior1(:,2,1,time) - normVals1(1))./normVals1(2);
            turnBiasInLightError = (flyError1(:,2,2,time))./normVals1(2);
            normVals2 = [nanmean(flyBehavior1(:,2,2,time)),nanstd(flyBehavior1(:,2,2,time))];
            turnBiasInDark = (flyBehavior1(:,2,2,time) - normVals2(1))./normVals2(2);
            turnBiasInDarkError = (flyError1(:,2,2,time))./normVals2(2);
            

            
            % Subtract expected standard error
            nullLDMErr = nanmean(turnBiasInLightError.^2 + turnBiasInDarkError.^2);
            
            % Calculate LDM by taking variance of LDM and subtracting null variance
            ldmDistribution = turnBiasInLight - turnBiasInDark;
            ldmCorrected_exp = sqrt(nanvar(ldmDistribution) - nullLDMErr);
            ldmCorrectedE_exp = bootstrp(numRe,@(x) sqrt(nanvar(x) - nullLDMErr),...
                ldmDistribution);
            
            % CONTROL GROUP
            
            % Set turn bias to variable
            normVals1 = [nanmean(flyBehavior2(:,2,1,time)),nanstd(flyBehavior2(:,2,1,time))];
            turnBiasInLight = (flyBehavior2(:,2,1,time) - normVals1(1))./normVals1(2);
            turnBiasInLightError = (flyError2(:,2,2,time))./normVals1(2);
            normVals2 = [nanmean(flyBehavior2(:,2,2,time)),nanstd(flyBehavior2(:,2,2,time))];
            turnBiasInDark = (flyBehavior2(:,2,2,time) - normVals2(1))./normVals2(2);
            turnBiasInDarkError = (flyError2(:,2,2,time))./normVals2(2);
      
            % Subtract expected standard error
            nullLDMErr = nanmean(turnBiasInLightError.^2 + turnBiasInDarkError.^2);
            
            % Calculate LDM by taking variance of LDM and subtracting null variance
            ldmDistribution = turnBiasInLight - turnBiasInDark;
            ldmCorrected_ctl = sqrt(nanvar(ldmDistribution) - nullLDMErr);
            ldmCorrectedE_ctl = bootstrp(numRe,@(x) sqrt(nanvar(x) - nullLDMErr),...
                ldmDistribution);
            
            
            ldmDeltaDist{c} = (ldmCorrectedE_exp - ldmCorrectedE_ctl)./ldmCorrectedE_ctl;
            ldm_experimental_groups{c} = [ldmCorrectedE_exp,ldmCorrectedE_ctl];

            
            geno_names{c} = ngn{ii};
            
            
            % Get N's for each experimental group
            n(c,:) = [sum(idx1),sum(idx2)]
            kk
            

            % Set group vector
            groups(c) = kk;
            c = c+1;

            
        
        
        end
    end



figure
line([0 length(ldmDeltaDist)+1],[0 0],'color',[0.5 0.5 0.5])

data = ldmDeltaDist;

for ii = 1:length(data)
        x = groups(ii)+rand*0.15;
        y_dat = median(data{ii});
        y_err = std(data{ii});
        hBar = errorbar(x,y_dat,y_err,'Marker','o','MarkerFaceColor','white','CapSize',0,...
            'MarkerSize',15,'LineWidth',1,'Color','black');
        hold on
end
hold off
shg



set(gca,'XTickLabel',names,'XTickLabelRotation',90)
xticks([1:23])
xlim([0 24])

ylabel('LDM')


pbaspect([5,3,1])
hFigure = gcf;
hFigure.Position =  [680   558   1600   600];
drawnow

medians = cellfun(@(x) nanmedian(x),data);

gal4s = geno_names;

if separate_by_driver
    
    eles = {'data_resamples','medians','neuron_names','n','groups','gal4s','eles'};
    eles_2 = {{'tb_light_experiment','tb_dark_experiment',...
        'turns_light_experiment','turns_dark_experiment'},{'tb_light_control','tb_dark_control',...
        'turns_light_control','turns_dark_control'},...
        eles{3:end}};
    if strcmp(eff1,'SHI')
        LDM_FigureS2C_Shi = {data,medians,names,n,groups,gal4s,eles};
        save('LDM_FigureS2C_Shi.mat','LDM_FigureS2C_Shi')       
        LDM_FigureS2AB_Shi= {bias_turns_experiment,bias_turns_control,eles_2};
        save('LDM_FigureS2AB_Shi.mat','LDM_FigureS2AB_Shi') 
    else
        LDM_FigureS2C_Trp = {data,medians,names,n,groups,gal4s,eles};
        save('LDM_FigureS2C_Trp.mat','LDM_FigureS2C_Trp')
        LDM_FigureS2AB_Trp= {bias_turns_experiment,bias_turns_control,eles_2};
        save('LDM_FigureS2AB_Trp.mat','LDM_FigureS2AB_Trp')
    end
    
else

    eles = {'data_resamples','medians','neuron_names','n','gal4s','eles'};
    if strcmp(eff1,'SHI')
        LDM_Figure2c = {data_temp,medians,names,n,gal4s,eles};
        save('LDM_Figure2c.mat','LDM_Figure2c')
    else
        LDM_Figure2d = {data,medians,names(groups),n,gal4s,eles};
        save('LDM_Figure2d.mat','LDM_Figure2d')
    end
    
end


%%

rh1 = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 20),'Un',0)'); 
rh7 = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 21),'Un',0)'); 
lpsp = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 12),'Un',0)'); 
epg = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 7),'Un',0)'); 
delta7 = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 14),'Un',0)');
pflc = cell2mat(cellfun(@(x) [median(x),std(x)],ldm_experimental_groups(groups == 11),'Un',0)');

raw = cell(0);
n_temp = cell(0);
for ii = 1:23
    raw{ii} = cell2mat(cellfun(@(x) [median(x),std(x)],...
        ldm_experimental_groups(groups == ii),'Un',0)'); 
    subplot(5,5,ii)
    dat = raw{ii};
    
    n_temp{ii} = n(groups==ii,:);
    
    for jj = 1:size(dat,1)
        errorbar(dat(jj,1:2),dat(jj,3:4))
        hold on
    end
    hold off
    xlim([0.5 2.5])
    
end

eles = {'data','error','names','groups','geno_names','n','eles'};

data = cellfun(@(x) x(:,1:2),raw,'Un',0);
error = cellfun(@(x) x(:,3:4),raw,'Un',0);

LDM_Figure2e = {data,error,names,groups,geno_names,n_temp,eles};
save('LDM_Figure2e.mat','LDM_Figure2e')


% subplot(1,6,1)
% dat = rh1;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% subplot(1,6,2)
% dat = rh7;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% subplot(1,6,3)
% dat = lpsp;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% 
% subplot(1,6,4)
% dat = epg;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% subplot(1,6,5)
% dat = delta7;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% 
% subplot(1,6,6)
% dat = pflc;
% for ii = 1:size(dat,1)
% errorbar(dat(ii,1:2),dat(ii,3:4))
% hold on
% end
% hold off
% xlim([0.5 2.5])
% 
% 
% shg

%% Figure 3d and 3e


% Sets neuron index
allLines = cell(0);
allBehave = cell(0);
        eff1 = 'SHI';
        eff2 = 'ISO';
        numRe = 10000;
         
actFilter = 0.02; % turns per second
flyN = nan(23,2);

allScores = cell(0);
scores_V1 = cell(0);
geno_names = cell(0);

n = nan(1,2);
groups = zeros(1);

separate_by_driver = true;
c = 1;
count = 1;

if separate_by_driver
    specific_lines = gal4s([28 40]);
else
    specific_lines = names([7 11 16 20]);
end
specific_line = cell(0);

for kk = 1:23
    %%
    
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
    
    if separate_by_driver
        num_drivers = length(ngn);
    else
        num_drivers = 1;
    end
    
    for ii = 1:num_drivers
        
        
        
        % Set index for specific experiment for shibire and control
        if separate_by_driver
            idx1 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff1);
            idx2 = ismember(upper(array(:,3)),ngn{ii}) & ismember(effectors,eff2);
            if sum(idx1) <= 50  || sum(idx2) <= 50
                continue
            end
            
        else
            idx1 = cell(0);
            idx2 = cell(0);
            for hh = 1:length(ngn)
                idx1{hh} = ismember(upper(array(:,3)),ngn{hh}) & ismember(effectors,eff1);
                idx2{hh} = ismember(upper(array(:,3)),ngn{hh}) & ismember(effectors,eff2);
                
                if sum(idx1{hh}) <= 50  || sum(idx2{hh}) <= 50
                    idx1{hh} = false(length(idx1{hh}),1);
                    idx2{hh} = idx1{hh};
                end
            end
            
            idx1 = any(cell2mat(idx1),2);
            idx2 = any(cell2mat(idx2),2);
            
            
            if sum(idx1) <= 50  || sum(idx2) <= 50
                continue
            end
        end
        
        
        
        ldmCorrected = nan(length(ngn),2);
        ldmCorrectedE = cell(0);
        ldmDeltaDist = cell(0);
        scores = cell(0);
        %         for ii = 1:length(ngn)
        
        
        
        % Extract turn biases for each genotype
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
        
        
        % V1
        scores = nan(numRe,4);
        for jj = 1:numRe
            index = randi(size(flyBiases1,1),size(flyBiases1,1),1);
            biases = flyBiases1(index,:);
            
            mean_shifts = mean(biases) - mean(biases)';
            %                     mean_shifts = zeros(size(mean_shifts));
            
            [a,b] = min(abs(biases(:,3) - biases(:,1:2) + mean_shifts(3,1:2)),...
                [],2);
            iLDM_exp_light = nanmean(b==1);
            [a,b] = min(abs(biases(:,4) - biases(:,1:2) + mean_shifts(4,1:2)),...
                [],2);
            iLDM_exp_dark = nanmean(b==2);
            
            
            index = randi(size(flyBiases2,1),size(flyBiases2,1),1);
            biases = flyBiases2(index,:);
            
            mean_shifts = mean(biases) - mean(biases)';
            %                     mean_shifts = zeros(size(mean_shifts));
            
            [a,b] = min(abs(biases(:,3) - biases(:,1:2) + mean_shifts(3,1:2)),...
                [],2);
            iLDM_ctl_light = nanmean(b==1);
            [a,b] = min(abs(biases(:,4) - biases(:,1:2) + mean_shifts(4,1:2)),...
                [],2);
            iLDM_ctl_dark = nanmean(b==2);
            
            scores(jj,:) = [iLDM_exp_light iLDM_exp_dark iLDM_ctl_light, iLDM_ctl_dark];
        end
        
        scores_V1{c} = scores;
        
        allScores{c} = scores;
        
        
        % V2
        scores = nan(numRe,8);
        for hh = 1:numRe
            indexes = randi(size(flyBiases1,1),size(flyBiases1,1),1);
            corrmap = corr(flyBiases1(indexes,:),'type','Spearman');
            fCorr = @(x,y) [x,y];
            corrs1 = [fCorr(corrmap(3,1),corrmap(3,2)), fCorr(corrmap(4,2),corrmap(4,1))];
            
            indexes = randi(size(flyBiases2,1),size(flyBiases2,1),1);
            corrmap = corr(flyBiases2(indexes,:),'type','Spearman');
            corrs2 = [fCorr(corrmap(3,1),corrmap(3,2)), fCorr(corrmap(4,2),corrmap(4,1))];
            
            scores(hh,:) = [corrs1,corrs2];
            
        end
        
        iLDM_light = sum(abs(corr(flyBiases1(:,1:2),flyBiases1(:,3),'type','Spearman')' - ...
            corr(flyBiases2(:,1:2),flyBiases2(:,3),'type','Spearman')'));
        iLDM_dark = sum(abs(corr(flyBiases1(:,1:2),flyBiases1(:,4),'type','Spearman')' - ...
            corr(flyBiases2(:,1:2),flyBiases2(:,4),'type','Spearman')'));
        
        scores = [iLDM_light iLDM_dark]
        scores_V2{c} = scores;
        kk
        
        % V3
        vectorMetric{c} = vectorBiasComponents({flyBiases1,flyBiases2});
        scores_V3(c) = vectorMetric(kk);
        
        
        geno_names{c} = ngn{ii};
        
        % Get N's for each experimental group
        n(c,:) = [sum(idx1),sum(idx2)]
        
        % Set group vector
        groups(c) = kk;
        c = c+1;
        
        
        
    end
    
    if ~separate_by_driver
        
        
        if kk == 7 || kk == 11 || kk == 20 || kk == 16
            specific_line{count} = {flyBiases1,flyBiases2};
            count = count+1;
        end
    else
        if kk == 40 || kk == 28
            specific_line{count} = {flyBiases1,flyBiases2};
            count = count+1;
        end
    end
    
end


clf

bar(cell2mat(cellfun(@(x) median(x),vectorMetric,'un',0)'))
xticklabels(names)
xticks(1:23)
shg



allScores_reshaped = nan(length(groups),4,numRe);
for ii = 1:length(allScores)
    try
        allScores_reshaped(ii,:,:) = permute(allScores{ii},[3,2,1]);
    end
end

geno_difference_scores = allScores_reshaped(:,[3 4],:) - allScores_reshaped(:,[1 2],:);



x_light = permute(geno_difference_scores(:,1,:),[3 1 2]);
x_dark =  permute(geno_difference_scores(:,2,:),[3 1 2]);

subplot(4,1,1)
try
    violin(x_light)
catch
    boxplot(x_light)
end
line([0 24],[0 0])

subplot(4,1,2)
try
    violin(x_dark)
catch
    boxplot(x_dark)
end
line([0 24],[0 0])


data_light = mat2cell(x_light,size(x_light,1),ones(size(x_light,2),1));
data_dark = mat2cell(x_dark,size(x_dark,1),ones(size(x_dark,2),1));

median_light = median(x_light);
median_dark = median(x_dark);


gal4s = geno_names;
subplot(4,1,3:4)
scatter(median_dark,median_light)
if separate_by_driver
    text(median_dark,median_light,gal4s)
else
text(median_dark,median_light,names(groups))
end
pbaspect([1 1 1])

eles = {'data_light_resamples','data_dark_resamples',...
    'medians_light','medians_dark','neuron_names','n','gal4s','eles'};

if separate_by_driver
    if strcmp(eff1,'SHI')
        LDM_FigureS3d_shi_drivers = {data_light,data_dark,median_light,median_dark,names,n,gal4s,eles};
        save('LDM_FigureS3d_shi_drivers.mat','LDM_FigureS3d_shi_drivers')
    elseif  strcmp(eff1,'TRP')
        LDM_FigureS3d_trp_drivers = {data_light,data_dark,median_light,median_dark,names(groups),n,gal4s,eles};
        save('LDM_FigureS3d_trp_drivers.mat','LDM_FigureS3d_trp_drivers')
    end
    
else
    if strcmp(eff1,'SHI')
        LDM_Figure3d_shi = {data_light,data_dark,median_light,median_dark,names,n,gal4s,eles};
        save('LDM_Figure3d_shi.mat','LDM_Figure3d_shi')
    elseif  strcmp(eff1,'TRP')
        LDM_Figure3d_trp = {data_light,data_dark,median_light,median_dark,names(groups),n,gal4s,eles};
        save('LDM_Figure3d_trp.mat','LDM_Figure3d_trp')
    end
end


%% Raw data for 3a

data = specific_line{4};
genotype = specific_lines(4);
eles = {'experiment_data','control_data','genotype','eles'};
LDM_Figure3a_rh1 = [data,genotype,eles];
save('LDM_Figure3a_rh1.mat','LDM_Figure3a_rh1')

data = specific_line{3};
genotype = specific_lines(3);
eles = {'experiment_data','control_data','genotype','eles'};
LDM_Figure3a_ibspsp = [data,genotype,eles];
save('LDM_Figure3a_ibspsp.mat','LDM_Figure3a_ibspsp')

data = specific_line{1};
genotype = specific_lines(1);
eles = {'experiment_data','control_data','genotype','eles'};
LDM_Figure3a_epg = [data,genotype,eles];
save('LDM_Figure3a_epg.mat','LDM_Figure3a_epg')

data = specific_line{2};
genotype = specific_lines(2);
eles = {'experiment_data','control_data','genotype','eles'};
LDM_Figure3a_pflc = [data,genotype,eles];
save('LDM_Figure3a_pflc.mat','LDM_Figure3a_pflc')

if separate_by_driver
   
    for ii = 1:2
        ngn = {'SS02252','SS00090'};
        idx1 = ismember(upper(array(:,3)),ngn(ii)) & ismember(effectors,eff1);
        idx2 = ismember(upper(array(:,3)),ngn(ii)) & ismember(effectors,eff2);
    
        % Extract turn biases for each genotype
        flyBehavior1 = cat(1,array{idx1,13});
        flyError1 = cat(1,array{idx1,14});
        activityI = flyBehavior1(:,1,2,1) > actFilter;
        flyBehavior1 = flyBehavior1(activityI,:,:,:);
        
        flyBehavior2 = cat(1,array{idx2,13});
        flyError2 = cat(1,array{idx2,14});
        activityI = flyBehavior2(:,1,2,1) > actFilter;
        flyBehavior2 = flyBehavior2(activityI,:,:,:);
        
        flyBiases1 = reshape(flyBehavior1(:,2,:,:),size(flyBehavior1,1),4); % Low Light, Low Dark, High Light, High Dark
        
        noNans = all(~isnan(flyBiases1),2);
        flyBiases1 = flyBiases1(noNans,:);
        
        flyBiases2 = reshape(flyBehavior2(:,2,:,:),size(flyBehavior2,1),4); % Low Light, Low Dark, High Light, High Dark
       
        noNans = all(~isnan(flyBiases2),2);
        flyBiases2 = flyBiases2(noNans,:);
        
        if ii == 1
            data = {flyBiases1,flyBiases2};
            genotype = specific_lines(2);
            eles = {'experiment_data','control_data','genotype','eles'};
            LDM_Figure3a_pflc = [data,genotype,eles];
            save('LDM_Figure3a_pflc.mat','LDM_Figure3a_pflc')
        else
            data = {flyBiases1,flyBiases2};
            genotype = specific_lines(1);
            eles = {'experiment_data','control_data','genotype','eles'};
            LDM_Figure3a_epg = [data,genotype,eles];
            save('LDM_Figure3a_epg.mat','LDM_Figure3a_epg')
        end
        
    end 
    
    
end


%%


for ii = 1:23
subplot(5,10,(ii-1)*2+1)
y_data = mean(allScores_reshaped(groups==ii,[1 3],:),3);
y_error = nanstd(allScores_reshaped(groups==ii,[1 3],:),[],3);
x_data = repmat([1,2],size(y_data,1),1);
x_data = x_data+normrnd(0,0.05,size(x_data));

color_1 = [0 0.5 1];
errorbar(x_data',y_data',y_error','Color',color_1,...
    'CapSize',0,'LineWidth',0.5,'Marker','o','MarkerSize',2,...
    'MarkerFaceColor',color_1,'MarkerEdgeColor','auto')
xlim([0 3])
ylim([0.4 0.8])
xticklabels({'Experiment','Control'})
xticks(1:2)
set(gca,'XTickLabelRotation',90)
title(names(ii))
shg

subplot(5,10,(ii-1)*2+2)
y_data = mean(allScores_reshaped(groups==ii,[2 4],:),3);
y_error = nanstd(allScores_reshaped(groups==ii,[2 4],:),[],3);
x_data = repmat([1,2],size(y_data,1),1);
x_data = x_data+normrnd(0,0.05,size(x_data));

color_1 = [1 0.15 0];
errorbar(x_data',y_data',y_error','Color',color_1,...
    'CapSize',0,'LineWidth',0.5,'Marker','o','MarkerSize',2,...
    'MarkerFaceColor',color_1,'MarkerEdgeColor','auto')
xlim([0 3])
ylim([0.4 0.8])
xticklabels({'Experiment','Control'})
xticks(1:2)
set(gca,'XTickLabelRotation',90)
title(names(ii))
shg
end
%%

eles = {'Light resamples','Dark Resamples','Light Medians','Dark Medians',"N's",'genotypes'};
outMat = {dataLight,dataDark,cellfun(@nanmedian,dataLight),cellfun(@nanmedian,dataDark),flyN,names,eles};
save('LDM_Figure3b-e','outMat')



