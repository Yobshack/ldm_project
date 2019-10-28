%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 05.24.2019 - updated
% Recreated for 60s on/off data for transitions


%% Select files to analyze

file_path = '../raw_data/60s_bins/';

% opens a ui to select mat files to analyze
fullPath = pwd;
files = uigetfile([file_path,'*.mat'],'MultiSelect','on');
if iscell(files) ~= 1
    files = {files};
end

tic
%turnDirs = batch(clust,@extractTurns,1,{files,lightDat});
toc

%turnDirs = cell(0);
tic
for ii = 1:length(files)
    f_name = [file_path,files{ii}];
    load(f_name) 
    turnDirs = extractTurns(flyTracks,f_name,lightDat);
    toc
    save([file_path,strcat(files{ii}(1:(end-4)),'_processed.mat')],'turnDirs')
ii
end
toc

%% Make large cell array to handle all data that I can query for groups

lightDat = dlmread('lightSequence_60s.txt');
lightDat = [lightDat;lightDat(end,:)];
lightDat(end,1) = lightDat(end,1)*10;

files = dir([file_path,'*processed.mat']);

[allDataCell, cellColNames] = parseProcessedFiles(files,file_path);


%% Generate a transition triggered average

array = allDataCell;

idx1 = strcmp(array(:,3),'CS');
idx2 = strcmp(array(:,3),'SOMA');
bin = 1;

idx = (idx1 | idx2);

idx = strcmp(array(:,3),'CS');

subArray = array(idx,:);

edgeMax = 60;
edge = 0:2:edgeMax;
tbs = nan(size(subArray,1),length(edge)-1,2);
err = tbs;
bias = nan(size(subArray,1),2,2);
minTurns = 10;
numRe = 100;

accum_time = cumsum(lightDat(:,1));
light = logical(lightDat(:,2));

ldm = nan(length(edge)-1,2,numRe);
for jj = 1:numRe
    
    for kk = 1:size(subArray,1)
        stamps = subArray{kk,2};
        index = randi(length(stamps),length(stamps),1);
        stamps = stamps(index);

        
        delta_times = accum_time'-stamps;
        delta_temp = cumsum(delta_times> 0,2);
        delta_temp = delta_temp == 1;
        
        
        turns = subArray{kk,1};
        turns = turns(index);
        
        idx1 = delta_temp & ~light';
        t_stamps{1} = delta_times(idx1);
        tr1 = turns(any(idx1,2));
        
        idx2 = delta_temp & light';
        t_stamps{2} = delta_times(idx2);
        tr2 = turns(any(idx2,2));
        
        [hc1,edge1,bin1] = histcounts(t_stamps{1},edge);
        
        [hc2,edge2,bin2] = histcounts(t_stamps{2},edge);
        
        for ii = 1:length(edge)-1
            if sum(bin1==ii) > minTurns && sum(bin2==ii) > minTurns
                
                tb = turnbias(tr1(bin1==ii,:));
                tbs(kk,ii,1) = tb(:,1);
                err(kk,ii,1) = tb(:,2);
                tb = turnbias(tr2(bin2==ii,:));
                tbs(kk,ii,2) = tb(:,1);
                err(kk,ii,2) = tb(:,2);
            end
        end
        
        bias(kk,:,1) = turnbias(turns(any(idx1,2)));
        bias(kk,:,2) = turnbias(turns(any(idx2,2)));
        
    end
    
    data = nanvar(tbs - bias(:,1,1),[],1);
    null = nanmean(err.^2 + bias(:,2,1).^2);
    
    ldm(:,:,jj) = squeeze(sqrt(data - null));
    jj
    
end

reshaped_ldm = abs(squeeze(cat(1,ldm(:,1,:),ldm(:,2,:))));
mean_ldm = nanmean(reshaped_ldm,2);

time = edge(1:(end-1))+diff(edge)/2;
time = [-fliplr(time) time]';
plot(time,mean_ldm,'-','Color',[1 0 0],'Marker','.')
patch_values = [mean_ldm+nanstd(reshaped_ldm,[],2);...
    flipud(mean_ldm-nanstd(reshaped_ldm,[],2))];
patch([time; flipud(time)],patch_values ,[1 0 0],'FaceAlpha',0.25,...
    'EdgeAlpha',0)
xlim([-15 15])
ylabel('LDM')
xlabel('Time Lag (s)')
shg

n = size(bias,1);
eles = {'dark_to_light_ldm','time_vector','n','eles'};
LDM_FigureS1F = {reshaped_ldm,time,n,eles};
save('LDM_FigureS1F.mat','LDM_FigureS1F')

n = size(bias,1);
eles = {'light_to_dark_ldm','time_vector','n','eles'};
LDM_FigureS1G = {reshaped_ldm([31:60 1:30],:),time,n,eles};
save('LDM_FigureS1G.mat','LDM_FigureS1G')

