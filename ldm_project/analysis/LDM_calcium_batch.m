%% batch ldm analysis file

startFolder = '/Volumes/LaCie/2p_experiments2/2*/2*';

% startFolder = '/Volumes/debivort_lab/Kyobi/Projects/2016.04.30_GCaMP_Imaging/raw_data/2p_experiments/2*/2*';
directories = dir(startFolder);


%% Analyze trials
%parpool(6)
failures = nan(1);
newIdx = 1;
for ii = 1:length(directories) 
    if directories(ii).isdir
        path = fullfile(directories(ii).folder,directories(ii).name);
        if isdir(fullfile(path,'tifs'))
            if ~isempty(ls(fullfile(path,'tifs')))
                try
                    disp(path)
                    cd(path)
                    load(strtrim(ls('*.mat')))
                    LDM_calcium_imaging(path,ballStruct)
                    ii
                catch
                    %             failures(newIdx) = ii;
                    %             newIdx = newIdx+1;
                    disp(ii)
                    
                end
            end
        end
    end 
    %c
end

%% Load individuals for further analysis
allData = cell(0);
avgInt = cell(0);
allPhys = allData;
c = 1;
clear indexes
index = 1:length(directories)
for ii = 1:length(directories)
    try
        indexes(c) = ii;
    path = fullfile(directories(index(ii)).folder,directories(index(ii)).name);
    cd(path)
    load(strtrim(ls('*state.mat')),'ballStruct');
    
    allData{c} = ballStruct.combined;
    allPhys{c} = ballStruct.physiology;
%     avgInt{c} = ballStruct.physiology.meanChannelValues;
    
    time = ballStruct.physiology.masterImageTime;
    
    c = c+1;
    
    catch
        
        ii
    
    end
    
commandLineProgressBar(ceil(ii/length(directories)),1,'!')
end

%%
gcamp = cell(30,1);
roll = gcamp;
pitch = gcamp;
yaws = gcamp;
corrmaps = gcamp;
for mm = 1:size(allData,2)
    mm
 
        yaws{mm} = allData{mm}.yaw;
%         if L_R_Index(mm)
%             yaws{mm} = -yaws{mm};
%         end
        pitch{mm} = allData{mm}.pitch;
        roll{mm} = allData{mm}.roll;
        gcamp{mm} = allData{mm}.GCaMP;
        corrmaps{mm} = allData{mm}.corrmap;
    
end


geno = {'2192'}; % change to 2252 to analyze PF-LC cells

subDirs = directories(indexes);
gIdx = contains({subDirs.name},geno) ;%& cellfun(@(x) ~isempty(x),allPhys);


ind = gIdx;
yaws = yaws(ind);
gcamp = gcamp(ind);
pitch = pitch(ind);

%% Output plot #1 - Overall response curve for many individuals


subAllPhys = allPhys(gIdx);
subAllData = allData(gIdx);

time = ballStruct.physiology.masterImageTime;


% This is the index of neurons that are right or left, 0 is left, 1 is right
L_R_Index = zeros(30,1);
L_R_Index([1, 12, 13, 14, 26]) = 1;
L_R_Index = logical(L_R_Index);
L_R_Index = L_R_Index(gIdx);



%%

pf_lc_index = contains({subDirs.name},'2252');
pf_lc_index = pf_lc_index(gIdx);

outCell = cell(0);

% Parameters

trials = 1:16

% Set indexes as conservative and to remove transition values
dIdx = time > -10 & time < 0;
dIdx1 = time > -50 & time < -10;
dIdx2 = time > 70 & time < 110;
lIdx = time > 10 & time < 50;

% number of resamples
numRe = 100;

mov_window = 25;
f = @(x) movmean((x - nanmean(x(:,dIdx),2))./nanmean(x(:,dIdx),2),mov_window,2);
% f = @(x) x./nanmean(x,2) - 1;

y1 = cell(length(subAllPhys),1);
y7 = y1;
y8 = y7;
y_corrmap = y1;
activity_all = y1;
clear sta staLD
for ii = 1:length(subAllPhys)
    try
       
        try
            y1{ii} = cat(1,subAllData{ii}.GCaMP(trials,:));
            y7{ii} = cat(1,subAllData{ii}.yaw(trials,:));
            y8{ii} = cat(1,subAllData{ii}.pitch(trials,:));
            y_corrmap{ii} = cat(3,subAllData{ii}.corrmap{trials});
        catch
            index = trials(1):...
                size(subAllData{ii}.GCaMP,1);
            y1{ii} = cat(1,subAllData{ii}.GCaMP(index,:));
            y7{ii} = cat(1,subAllData{ii}.yaw(index,:));
            y8{ii} = cat(1,subAllData{ii}.pitch(index,:));
            
            try
            y_corrmap{ii} = cat(3,subAllData{ii}.corrmap{index});
            end
        end
       
        
        
        activity = nansum(abs(y7{ii}),2);
        idx = activity<250 & activity > 15*pi;
        activity_all{ii} = activity;
        
        fluor = nanmean(y1{ii},2);
        minFluor = 10;
        idx = idx & fluor > minFluor;
        
        idx = idx & nanstd(y1{ii},[],2)./nanmean(y1{ii},2) < 0.8;
        
        y1{ii} = y1{ii}(idx,:);
        y7{ii} = y7{ii}(idx,:);
        y8{ii} = y8{ii}(idx,:);
        try
        y_corrmap{ii} = y_corrmap{ii}(:,:,idx);
        end
        activity_all{ii} = activity_all{ii}(idx);
        
        clf
        colors = jet(size(y1{ii},1));
        hPlot = plot(time,y1{ii});
        for kk = 1:length(hPlot)
            hPlot(kk).Color = [colors(kk,:) 0.5];
        end
        hold on
        plot(time,nanmean(y1{ii}),'LineWidth',3,'Color',[0 0 0])
        %     y2{ii} = nanmedian(y1{ii});
        shg
        %     pause(2)
        
%         sta(ii,:,:) = cat(3,subAllData{ii}.rightTurnSTA,subAllData{ii}.leftTurnSTA);
%         staLD(ii,:,:) = cat(3,subAllData{ii}.allTurnSTALight,subAllData{ii}.allTurnSTADark);
%         turnVals{ii} = subAllData{ii}.turnVals;
    end
end
shg

%%

SR = 1/(time(2)-time(1))
min_trials = 1;
idx1 = cellfun(@(x) size(x,1) >= min_trials,y1);
idx2 = cellfun(@(x) size(x,1) >= min_trials,y7);
y3 = cellfun(@(x) f(x),y1(idx1),'Un',0);
y1 = y1(idx1);
y7 = cellfun(@(x) movmean(x/SR,mov_window,2),y7(idx2),'Un',0);
y8 = cellfun(@(x) movmean(x/SR,mov_window,2),y8(idx2),'Un',0);
activity_all = activity_all(idx2);
% y7 = cellfun(@(x) movmean(x,20,2),y7,'un',0);
y4 = cell2mat(cellfun(@(x) nanmedian(x,1),y3,'Un',0));
y5 = nanstd(y4,[],1);
y6 = cell2mat(cellfun(@(x) nanmean(x,1),y1,'Un',0));

y_corrmaps = y_corrmap(idx1);


L_R_Index_sub = L_R_Index(idx1);

pf_lc_index = pf_lc_index(idx1);


%% Make matrices for Ben

yaw = cell2mat(y7);
pitch = cell2mat(y8);
gcamp = cell2mat(y3);

clear ind_index
for ii = 1:size(y7,1)
    ind_index{ii} = repmat(ii,size(y7{ii},1),1);

end

ind_index = cell2mat(ind_index');

eles = 'time,gcamp,yaw,pitch,index,eles';
LDM_Figure4GH = {time,gcamp,yaw,pitch,ind_index,eles,geno};
save('LDM_Figure4GH.mat','LDM_Figure4GH')


%% Plots for paper
clf
% y7 = cellfun(@(x) x.yaw,subAllData(idx1),'un',0);

f_1 = @(x) mean(x,2);

f_2 = @(x) sum(sum(abs(x(x<=0)),2))./sum(sum(abs(x),2),1); % Makes proportion right turn

% f_2 = sum(x.*(x<0.0001),2)./sum(abs(x.*(abs(x)>0.0001)),2)

% f_2 = @(x) sum(cumsum(x.*(x<0),2),2)./sum(cumsum(abs(x),2),2); 

d1 = cell2mat(cellfun(@(x) [f_1(x(:,dIdx1)),...
    f_1(x(:,lIdx)),f_1(x(:,dIdx2))],y3,'un',0));

d2 = cell2mat(cellfun(@(x) [f_2(x(:,dIdx1)),...
    f_2(x(:,lIdx)),f_2(x(:,dIdx2))],y7,'un',0));


numRe = 10000;

corrs1 = bootstrp(numRe,@(x,y) corr(x,y,'type','Spearman'),d2(:,1),d2(:,3));
corrs2 = bootstrp(numRe,@(x,y) corr(x,y,'type','Spearman'),d2(:,1),d2(:,2));


subplot(2,2,1)
if size(d2,1) > 20
n = 23
else
    n = 1;
end
plot(d2(n,:),'color','black')
xlim([0.5 3.5])
ylim([0 1])
ylabel('Proportion right rotation')


subplot(2,2,2)
scatter(d2(:,1),d2(:,3),'Marker','o','MarkerFaceColor','flat','CData',double(pf_lc_index));

% colormap(winter)
% ylim([0 1])
% xlim([0 1])

subplot(2,2,3)
hold on
scatter(d2(:,1),d2(:,2),'Marker','o','MarkerFaceColor','flat','CData',double(pf_lc_index));
% 
% ylim([0 1])
% xlim([0 1])

subplot(2,2,4)
violin([corrs1,corrs2])
xticklabels({'DD','LD'})
xticks(1:2)

shg

if size(d2,1) >= 20
    genos = cellstr(repmat('2252',size(d1,1),1));
    genos(~pf_lc_index) = {'2192'};
    individual = d2(23,:);
    scatter_1 = [d2(:,1),d2(:,3)];
    scatter_2 = [d2(:,1),d2(:,2)];
    eles = {'individual','scatter_DD','scatter_LD','correlation_resamples','n','genos','eles'}
    LDM_Figure4b = {individual,scatter_1,scatter_2,[corrs1,corrs2],size(d1,1),genos,eles};
    save('LDM_Figure4b.mat','LDM_Figure4b');
end


