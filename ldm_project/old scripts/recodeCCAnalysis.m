%% Select files to analyze

% opens a ui to select mat files to analyze
fullPath = pwd;
files = uigetfile('*.mat','MultiSelect','on');
if iscell(files) ~= 1
    files = {files};
end

% loads an experimental design file, script expects a n x 3 array where n is the number of bins in
% the stimulus sequence. col 1 is the bin length in seconds, col 2 is a logical vector indicating
% whether the lights were on or off in that bin. col 3 is a vector of pwm values that reflect the
% pwm setting of the arduino controlling the led light board.
lightDat = dlmread('lightSequenceScreen.txt');
lightDat = [lightDat;lightDat(end,:)];
lightDat(end,1) = lightDat(end,1)*10;
tic

%turnDirs = cell(0);
tic
for ii = 1:length(files)
    load(files{ii}) 
    turnDirs{ii} = extractTurns(flyTracks,files{ii},lightDat);
    toc
    save(strcat(files{ii}(1:(end-4)),'_processed.mat'),'turnDirs')
end
toc


%% Analyze LDM from all experiments

tfilter = 20;
numRe = 1;
timeLow = [0 60]*60;
timeHigh = [120 180]*60;

rubinAnno = importsplit('/Users/kyobikakaria/Desktop/Data/split_screen/split_Rubin_Annotation.csv');
neuronNames = table2array(rubinAnno(:,1));
splitNames = rubinAnno.Properties.VariableNames(2:end);
tic

ldmRaw = cell(0); lightBias = ldmRaw; lightTurns = ldmRaw; darkBias = ldmRaw; darkTurns = ldmRaw;
ldm = ldmRaw; ldmSim = ldmRaw; correctedLDM = ldmRaw; correctedLDMErr = ldmRaw;
intraChange = ldmRaw; intraChangeErr = ldmRaw; allLightLike = ldmRaw;
allLightLikeE = ldmRaw; allDarkLike = ldmRaw;
allDarkLikeE = ldmRaw; switchDarkMat = ldmRaw; switchDarkMatE = ldmRaw; switchLightMat = ldmRaw;
switchLightMatE = ldmRaw;

clear sl sd sle sde data2 error2 sl2 sd2 sle2 sde2 geno effector neuron correctedLDMRe lightBar lightBarE
clear darkBar darkBarE



for jj = 1:length(turnDirs)
    jj
    label = strsplit(char(turnDirs{jj}.tLabels(1,1)),'_');
    geno(jj) = label(1);
    effector(jj) = upper(label(2));
    nName = neuronNames(logical(sum(table2array(rubinAnno(:,strmatch(geno{jj},splitNames)+1)),2)),:);
    if ~isempty(nName)
        neuron(jj) = {[nName{:}]};
    else
        neuron(jj) = {''};
    end
    
    for ii = 1:size(turnDirs{jj}.tSequence)
        fly = ii;
        
        
        
        turns = turnDirs{jj}.tSequence{fly};
        times = turnDirs{jj}.tTime{fly};
        
        cumTime = [0;cumsum(lightDat(:,1))];
        
        lights = ones(length(times),1);
        for kk = 1:length(times)
            
            lights(kk,1) = lightDat(sum(cumTime < times(kk)),2);
        end
        
        lightBias{jj}(fly,1) = nanmean(turns(times>timeLow(1) & times<=timeLow(2) & lights));
        lightBias{jj}(fly,2) = nanmean(turns(times>timeHigh(1) & times <=timeHigh(2) & lights));
        
        lightTurns{jj}(fly,1) = length(turns(times>timeLow(1) & times<=timeLow(2) & lights));
        lightTurns{jj}(fly,2) = length(turns(times>timeHigh(1) & times <=timeHigh(2) & lights));
        
        darkBias{jj}(fly,1) = nanmean(turns(times>timeLow(1) & times <=timeLow(2) & ~lights));
        darkBias{jj}(fly,2) = nanmean(turns(times>timeHigh(1) & times <=timeHigh(2) & ~lights));
        
        darkTurns{jj}(fly,1) = length(turns(times>timeLow(1) & times <=timeLow(2) & ~lights));
        darkTurns{jj}(fly,2) = length(turns(times>timeHigh(1) & times <=timeHigh(2) & ~lights));
        
        genoName{jj}(fly) = geno(jj);
        effectorName{jj}(fly) = effector(jj);
        neuronName{jj}(fly) = neuron(jj);
    end
    
    
    lightBias{jj}(lightTurns{jj} <= tfilter) = nan;
    darkBias{jj}(lightTurns{jj} <= tfilter) = nan;
    
    ldmRaw{jj} = lightBias{jj} - darkBias{jj};
    ldm{jj} = nanmean(abs(ldmRaw{jj}));
    
    avgBias = nanmean(cat(3,lightBias{jj},darkBias{jj}),3);
    
    ldmSim{jj} = nanmean(nanmean(abs((binornd(repmat(lightTurns{jj},1,1,500),repmat(avgBias,1,1,500))./...
        repmat(lightTurns{jj},1,1,500)) - (binornd(repmat(darkTurns{jj},1,1,500),...
        repmat(avgBias,1,1,500))./repmat(darkTurns{jj},1,1,500))),1),3);
    
    correctedLDM{jj} = ldm{jj} - ldmSim{jj};
    
    for kk = 1:numRe
        idx = randi(size(ldmRaw{jj},1),size(ldmRaw{jj},1),1);
        ldmRe = ldmRaw{jj}(idx,:);
        ldmSimRe = nanmean(nanmean(abs((binornd(repmat(lightTurns{jj}(idx,:),1,1,10),...
            repmat(avgBias(idx,:),1,1,10))./...
            repmat(lightTurns{jj}(idx,:),1,1,10)) - (binornd(repmat(darkTurns{jj}(idx,:),1,1,10),...
            repmat(avgBias(idx,:),1,1,10))./repmat(darkTurns{jj}(idx,:),1,1,10))),1),3);
        correctedLDMRe(kk,:) = nanmean(abs(ldmRe(idx,:))) - ldmSimRe;
    end
    
    correctedLDMErr{jj} = nanstd(correctedLDMRe,1);
    
    intraChange{jj} = diff(correctedLDM{jj})./nanmean(abs(ldm{jj}(:,1)));
    intraChangeErr{jj} = sqrt(sum(correctedLDMErr{jj}.^2))/(length(correctedLDMErr{jj})*...
        nanmean(abs(ldm{jj}(:,1))));
    
    
    % Calculate light/dark shifts
    lightTBs = lightBias{jj};
    darkTBs = darkBias{jj};
    
    if all(isnan(lightTBs(:,1)))
        continue
    end
    
    x1 = lightTBs(:,2); x2 = darkTBs(:,2); x3 = lightTBs(:,1); x4 = darkTBs(:,1);
    
    [a,b] = min(abs([(x1 - x3) (x1 - x4)]),[],2);
    
    lightBar = nansum(b(~isnan(a))==1)/length(b(~isnan(a))) - 0.5;
    
    [a,b] = min(abs([(x2 - x3) (x2 - x4)]),[],2);
    darkBar =  nansum(b(~isnan(a))==1)/length(b(~isnan(a))) - 0.5;
    
    
    for ii = 1:100
        rind = randi(length(x1),length(x1),1);
        
        [a,b] = min(abs([(x1(rind) - x3(rind)) (x1(rind) - x4(rind))]),[],2);
        lightBarE(ii) = nansum(b(~isnan(a))==1)/length(b(~isnan(a))) - 0.5;
        
        [a,b] = min(abs([(x2(rind) - x3(rind)) (x2(rind) - x4(rind))]),[],2);
        darkBarE(ii) =  nansum(b(~isnan(a))==1)/length(b(~isnan(a))) - 0.5;
    end
    
    
    allLightLike{jj} = lightBar;
    allDarkLike{jj} = darkBar;
    allLightLikeE{jj} = std(lightBarE);
    allDarkLikeE{jj} = std(darkBarE);
    
end
toc
%% Group by neurons
correctedMat = cell(length(neuronNames),3);
correctedMatE = correctedMat;
intraMat = correctedMat;
intraMatE = correctedMat;
effectorTypes = unique(effector);
offsets = [-0.25 0 0.25];
hold on

color = colormap(gray(16));
for ii = 1:length(neuronNames)
    clear data error
    for jj = 1:2
        index = strcmp(neuron,neuronNames(ii)) &...
            strcmp(effector,effectorTypes(jj));
        splitName = geno(index);
        [a,b,c] = unique(splitName);
                correctedMat{ii,jj} = correctedLDM(strcmp(neuron,neuronNames(ii)) &...
                    strcmp(effector,effectorTypes(jj)));
                correctedMatE{ii,jj} = correctedLDMErr(strcmp(neuronNames(ii),neuron) &...
                    strcmp(effector,effectorTypes(jj)));
%                 if ~isempty(correctedMat{ii,jj})
%                     data = cat(1,correctedMat{ii,jj}{:});
%                     error = cat(1,correctedMatE{ii,jj}{:});
%                     bar(ii+offsets(jj),nanmean(data(:,2)),'BarWidth',0.1)
%                     errorbar(ii+offsets(jj),nanmean(data(:,2)),sqrt(sum(error(:,2).^2))/size(data,1))
%                     for kk = 1:length(correctedMat{ii,jj})
%                         plot(ii+offsets(jj),correctedMat{ii,jj}{kk}(:,2),'Marker','.','MarkerSize',30,...
%                             'Color',color(kk,:))
%                     end
%                 end
        
        intraMat{ii,jj} = intraChange(index);
        intraMatE{ii,jj} = intraChangeErr(index);
        
        switchLightMat{ii,jj} = allLightLike(index);
        switchDarkMat{ii,jj} = allDarkLike(index);
        switchLightMatE{ii,jj} = allLightLikeE(index);
        switchDarkMatE{ii,jj} = allDarkLikeE(index);
        
        if ~isempty(intraMat{ii,jj})
            data{jj} = cat(1,intraMat{ii,jj}{:});
            error{jj} = cat(1,intraMatE{ii,jj}{:});
            
            sl{jj} = cat(1,switchLightMat{ii,jj}{:});
            sd{jj} = cat(1,switchDarkMat{ii,jj}{:});
            
            sle{jj} = cat(1,switchLightMatE{ii,jj}{:});
            sde{jj} = cat(1,switchDarkMatE{ii,jj}{:});
            
            if jj > 1
                data2(ii) = (nanmean(data{jj}) - nanmean(data{1}));
                error2(ii) = sqrt(nanmean(error{jj}.^2)/length(error{jj}).^2 +...
                    nanmean(error{1}.^2)/length(error{1}.^2));
                
                sl2(ii) = (nanmean(sl{jj}) - nanmean(sl{1}));
                sle2(ii) = sqrt(nanmean(sle{jj}.^2)/length(sle{jj}).^2 +...
                    nanmean(sle{1}.^2)/length(sle{1}.^2));
                
                sd2(ii) = (nanmean(sd{jj}) - nanmean(sd{1}));
                sde2(ii) = sqrt(nanmean(sde{jj}.^2)/length(sde{jj}).^2 +...
                    nanmean(sde{1}.^2)/length(sde{1}.^2));
                
                
                
                %             for kk = 1:length(a)
                %                 plot(ii+offsets(jj),(nanmean(data{jj}(c==kk)) - ...
                %                     nanmean(data{1}(c==kk)))*-1,...
                %                     'Marker','.','MarkerSize',30,...
                %                     'Color',color(kk,:))
                %             end
                
                
            end
        end
    end
end

%%

hold on

labs = {'P-FN1','P-FN2','P-FN3','P-FN4','P-EN','P-EG','E-PG','E-PG_Cap','PF-FG','PF-FRub',...
    'PF-LCre','LPs-P','PBint-Cap','PBint-d7','Ps-P','PsIb-P'};
b1 = bar(data2,'BarWidth',0.5,'FaceColor',[1 0 0]);
e = errorbar(data2,error2,'Color',[0 0 0],'LineWidth',3,'LineStyle','none');

set(gca,'XTickLabel',labs,'XTickLabelRotation',45,'XTick',1:length(neuronNames),...
    'FontWeight','bold','FontSize',16)


for ii = 1:length(data2)
    d1 = data2(ii);
    e1 = error2(ii);
    if isnan(d1)
        continue
    end
    pd = makedist('Normal',abs(d1),e1);
    pVal = cdf(pd,0);
    if pVal > 0.05
        character = '';
    elseif pVal <= 0.05 && pVal > 0.01
        character = '*';
    elseif pVal <= 0.01 && pVal > 0.001
        character = '**';
    elseif pVal <= 0.001
        character = '***';
    end
    text(b1(jj).XData(ii) + b1(jj).XOffset,...
        b1(jj).YData(ii)+(e(jj).YPositiveDelta(ii)*1.6*(b1(jj).YData(ii)/abs(b1(jj).YData(ii)))),...
        character,'FontSize',16,'HorizontalAlignment','center')
    drawnow
end

shg

%% Calculating "light-like" shifts
clf
hold on

colormap(parula)
labs = {'P-FN1','P-FN2','P-FN3','P-FN4','P-EN','P-EG','E-PG','E-PG_Cap','PF-FG','PF-FRub',...
    'PF-LCre','LPs-P','PBint-Cap','PBint-d7','Ps-P','PsIb-P'};
b1 = bar([sl2 ;sd2]','BarWidth',0.75);
drawnow
e = errorbar([(1:length(sl2))+b1(1).XOffset ; (1:length(sl2))+b1(2).XOffset]',...
    [sl2 ;sd2]',[sde2 ;sde2]','Color',[0 0 0],'LineWidth',2,'LineStyle','none','CapSize',0);

set(gca,'XTickLabel',labs,'XTickLabelRotation',45,'XTick',1:length(neuronNames),...
    'FontWeight','bold','FontSize',16)

data3 = [sl2',sd2'];
error3 = [sle2',sde2'];

for ii = 1:length(data2)
    for jj = 1:2
    d1 = data3(ii,jj);
    e1 = error3(ii,jj);
    if isnan(d1)
        continue
    end
    pd = makedist('Normal',abs(d1),e1);
    pVal = cdf(pd,0);
    if pVal > 0.05
        character = '';
    elseif pVal <= 0.05 && pVal > 0.01
        character = '*';
    elseif pVal <= 0.01 && pVal > 0.001
        character = '**';
    elseif pVal <= 0.001
        character = '***';
    end
    text(b1(jj).XData(ii) + b1(jj).XOffset,...
        b1(jj).YData(ii)+(e(jj).YPositiveDelta(ii)*1.6*(b1(jj).YData(ii)/abs(b1(jj).YData(ii)))),...
        character,'FontSize',16,'HorizontalAlignment','center')
    drawnow
    end
end

legend('light','dark')

shg
