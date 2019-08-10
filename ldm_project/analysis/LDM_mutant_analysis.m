%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 05.23.2019 - updated


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
lightDat = dlmread('lightSequence_5Mins.txt');
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
[allDataCell, cellColNames] = parseProcessedFiles(files);

%%
    
array = allDataCell;

array(:,3) = upper(array(:,3));


%% Generate ldm bar plot

[a,b,c] = unique(array(:,3));

ldmNull = nan(1,2);
numRe = 10000;
ldm = nan(length(a),numRe);
f = @(x) (x - nanmean(x))./nanstd(x);
for ii = 1:length(a)
    turns = array(c==ii,1);
    light = array(c==ii,11);
   
    dat = nan(sum(c==ii),1);
    err = dat;
    for jj = 1:length(turns);
        tb_light = turnbias(turns{jj}(light{jj}));
        tb_dark = turnbias(turns{jj}(~light{jj}));
        dat(jj) = tb_light(:,1) - tb_dark(:,1);
        err(jj) = sqrt(tb_light(:,2).^2 + tb_dark(:,2).^2);
    end
    index = ~isnan(dat);
    dat = dat(index);
    err = err(index);

    ldm(ii,:) = abs(bootstrp(numRe,@(dat,err) sqrt(var(dat) - mean(err.^2)),dat,err)');
    n(ii) = min(sum(~isnan(dat)));
end

id = strcmp(a(e),'ISO31');
[d,e] = sort(nanmean(ldm,2));
ldm = ldm(e,:);
n = n(e);

names = {'CS','cry[1]','cry[2]','Rh4','eya',...
    'gmr-hid','iso31','Rh1[7]','Rh1[8]','NorpA-'};

errorbar(nanmean(ldm,2),nanstd(ldm,[],2),'LineStyle','none','CapSize',0,'MarkerSize',10,'Marker','o',...
    'MarkerFaceColor','auto','LineWidth',2)
set(gca,'XLim',[0 11],'YLim',[0 0.15],'XTickLabel',names(e),'XTickLabelRotation',45,...
    'XTick',1:10)
ylabel('LDM')
pbaspect([1.5 1 1])
eles = {'ldm_resamples','names','n','eles'};
LDM_Figure_2a = {ldm,names(e),n(e),eles};
save('LDM_Figure_2a.mat','LDM_Figure_2a');


shg
