%% Kyobi Skutt-Kakaria 
% Bee and Cricket Analysis
% 11.20.2018 - Created
% 04.27.2019 - Updated


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
% lightDat = dlmread('lightSequenceScreen_15Mins.txt');
% lightDat = [lightDat;lightDat(end,:)];
% lightDat(end,1) = lightDat(end,1)*10;
tic

lightDat = [repmat(900,96,1), repmat([1;0],48,1),repmat([75;0],48,1)];
%turnDirs = batch(clust,@extractTurns,1,{files,lightDat});
toc

flyTracksCell = cell(0);
tic
for ii = 1:length(files)
    rootName = files{ii}(1:end-4);
    load(files{ii})
    
    timeT = fopen(strcat(rootName,'_Time.bin'));
    flyTracks.tStamps = cumsum(fread(timeT,'double'));   
    
    turnT = fopen(strcat(rootName,'_Turns.bin'));
    flyTracks.rightTurns = fread(turnT,[10 length(flyTracks.tStamps)],'double');   
    
    flyTracks.nFlies = size(expmt.labels_table,1);
    
    flyTracks.mazeOri = [zeros(1,6),ones(1,4)];
    
    flyTracks.lightDat = lightDat;

    turnDirs = calculateTurnDirsV2(flyTracks);
    
    flyTracks.turnDirs = turnDirs;
    
    flyTracksCell{ii} = flyTracks;
    
end
toc

        accumTime = cumsum(flyTracksCell{ii}.lightDat(:,1));
%% Analyze LDM for bees and crickets

tb = cell(cell(0));
turnC = tb;
for ii = 1:length(files)
   ii
    

    clf
    for jj = 1:10
        
        % defines how many bins the analysis will have
        binSum = 1 * 60*60;
        
        
        % returns time vector from calculate turn dirs
        time = flyTracksCell{ii}.turnDirs.tTime{jj};
        if length(time) < 100
            continue
        end
        
        
        % Defines what time bin to consider
        totalTimeCutoff = [0 24]*60*60; % time in hours converted to seconds
        idx = time > totalTimeCutoff(1) & time < totalTimeCutoff(2);
        time = time(idx);
        
        % Taking the mod of the time by the binSum combines the data from all hours
        tstamps = mod(time,binSum);  
        binLength = 15*60;
        turns = flyTracksCell{ii}.turnDirs.tSequence{jj};
        turns = turns(idx,:);

        turnN = length(turns);
        
 
        
        for kk = 1:round(binSum/binLength)
            

            idx = tstamps>(kk-1) * binLength & tstamps < (kk) * binLength;
            tb{ii}{jj}(kk,:) = turnbias(turns(idx));
            turnC{ii}{jj}(kk,:) = length(turns(idx));

        end
        
        
        
        subplot(5,2,jj)
        errorbar(tb{ii}{jj}(:,1),tb{ii}{jj}(:,2),'CapSize',0,'Marker','.','MarkerSize',5)
        
        ylim([0.25 0.75])
        xlim([0 size(tb{ii}{jj},1)])
        
        line([0 9],[0.5 0.5],'color','red')
        
        
    end
    shg
    pause(1)
end

%%
mat = zeros(1,4,2);
animal = zeros(1);
c = 1
for ii = 1:size(tb,2)
    ii
    
    for jj = 1:size(tb{ii},2)
        try
            mat(c,:,:) = reshape(tb{ii}{jj},4,1,2);
            matT(c,:) = turnC{ii}{jj}';
            
            if ii == 1 || ii == 2
                animal(c) = 1; % cricket
            elseif ii == 3 || ii == 4
                animal(c) = 2; % bee
            end
            c = c+1;
        end
    end
    
    
end

idx = all(matT>20,2);
mat = mat(idx,:,:)
animal = animal(idx);


colors = [0, 1, 0;
    0, 0, 1]
colormap(colors)

subplot(1,4,1)
hScat = scatter(mat(:,1,1),mat(:,3,1),'CData',animal);
hScat.MarkerFaceColor = hScat.MarkerEdgeColor;
xlim([0.25 0.75])
ylim([0.25 0.75])
xlabel('Turn Bias in Light')
ylabel('Turn Bias in Light')
pbaspect([1 1 1])

subplot(1,4,2)
hScat = scatter(mat(:,2,1),mat(:,4,1),'CData',animal);
hScat.MarkerFaceColor = hScat.MarkerEdgeColor;
xlim([0.25 0.75])
ylim([0.25 0.75])
xlabel('Turn Bias in Dark')
ylabel('Turn Bias in Dark')
pbaspect([1 1 1])


subplot(1,4,3)
avgMat = [mean(mat(:,[1 3],1),2), mean(mat(:,[2 4],1),2)];
hScat = scatter(avgMat(:,1),avgMat(:,2),'CData',animal);
hScat.MarkerFaceColor = hScat.MarkerEdgeColor;
xlim([0.25 0.75])
ylim([0.25 0.75])
xlabel('Turn Bias in Dark')
ylabel('Turn Bias in Light')
pbaspect([1 1 1])

subplot(1,4,4)
numRe = 100;
corrs = corr(mat(:,:,1),'type','Spearman');
corrs = [corrs(1,3),corrs(2,4),corrs(1,2)];
corrsE = nan(numRe,size(corrs,2));
for ii = 1:numRe
    ind = randi(size(mat,1),size(mat,1),1);
    c = corr(mat(ind,:,1),'type','Spearman');
    corrsE(ii,:) = [c(1,3),c(2,4),c(1,2)];
end
hBar = bar(corrs);
hBar.FaceColor = [0.5 0.5 0.5];
set(gca,'XTickLabel',{'LvL','DvD','LvD'})
ylabel('Spearman Correlation Coefficient')
hold on
hError = errorbar(hBar.XData,hBar.YData, std(corrsE,[],1),'LineStyle','none','CapSize',0,...
    'LineWidth',3);
hold off
hError.Color = [0 0 0];
pbaspect([1 1 1])


[bestAn,I] = sortrows([avgMat(:,1) - avgMat(:,2),animal']);

figure
time = [15 30 45 60];
one = I(23);
two = I(24);
hPlot = plot(time,[mat(one,:,1);mat(two,:,1)]','LineWidth',3);
hold on
hError = errorbar(repmat(time,2,1)',[mat(one,:,1);mat(two,:,1)]',[mat(one,:,2);mat(two,:,2)]','LineStyle','none','CapSize',0,...
    'LineWidth',3,'Color',[0 0 0]);
hold off
for ii = 1:length(hPlot)
    hPlot(ii).Color = colors(ii,:);
    hError(ii).Color  = colors(ii,:);
end

legend('Cricket','BumbleBee')
legend('boxoff')
ylim([0 1])
xlim([0 75])
pbaspect([2,1,1])
ylabel('Turn Bias')
xlabel('Time in minutes')

n = size(mat,1);
eles = {'light_light','dark_dark','light_dark','resampled_correlation','bumble_cricket_index',...
    'individual_biases','individual_binomial_standarderror','n'};
LDM_Figure1j = {[mat(:,1,1),...
    mat(:,3,1)],[mat(:,2,1),mat(:,4,1)],...
    avgMat,corrsE,animal,[mat(one,:,1);mat(two,:,1)]',...
    [mat(one,:,2);mat(two,:,2)]',n,eles};

save('LDM_Figure1j.mat','LDM_Figure1j')

