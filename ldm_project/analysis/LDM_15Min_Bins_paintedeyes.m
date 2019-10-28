%% Kyobi Skutt-Kakaria
% 05.02.2019 - Created


%% Select files to analyze

% opens a ui to select mat files to analyze
fullPath = pwd;
file_path = '~/git/ldm_project/raw_data/painted_eyes/';
files = dir([file_path,'*.mat']);

lightDat = dlmread('~/git/ldm_project/resources/lightSequenceScreen_15Mins.txt');
lightDat = [lightDat;lightDat(end,:)];
lightDat(end,1) = lightDat(end,1)*10;

%turnDirs = cell(0);
for ii = 1:length(files)
    f_name_raw = files(ii).name;
    f_name = [file_path,f_name_raw];
    load(f_name) 
    turnDirs = extractTurns(flyTracks,f_name,lightDat);
    toc
    save([file_path,strcat(f_name_raw(1:(end-4)),'_processed.mat')],'turnDirs')

end
toc


%% Make large cell array to handle all data that I can query for groups

files = dir([file_path,'*processed.mat']);
[allDataCell, cellColNames] = parseProcessedFiles(files,file_path);



%%
array = allDataCell;

array(:,3) = upper(array(:,3));


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

%%
hErr = errorbar(1:3,mean(ldm,2),std(ldm,[],2),'LineStyle','none','Marker','s',...
'color','k');
hErr.MarkerFaceColor = hErr.MarkerEdgeColor;
hErr.MarkerSize = 10;
hErr.CapSize = 7.5;
hErr.LineWidth = 1.5;
xlim([0.5 3.5])
ylim([0 0.15])
text(hErr.XData,hErr.YData - 0.05,num2str(n))
shg
