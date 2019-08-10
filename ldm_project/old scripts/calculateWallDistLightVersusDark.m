lightDat = dlmread('lightSequenceScreen.txt');
lightDat = [lightDat;lightDat(end,:)];
lightDat(end,1) = lightDat(end,1)*10;

%%
allFliesDist = cell(length(files),1);
for ii = 2:4
    ii
    load(files{ii})
    flyTracks = adjustTimeStamps(flyTracks);
    flyTracks = distFromWall(flyTracks,lightDat);
    allFliesDist{ii} = [flyTracks.distLight flyTracks.distDark];
end

allDist = cat(1,allFliesDist{:});
[~,p] = ttest(allDist(:,1),allDist(:,2));
%% Make figure of distance to wall
clf
histogram(allDist(:,2)+2,'Normalization','probability','FaceColor',[0 0 0.5])
hold on
histogram(allDist(:,1)+2,'Normalization','probability','FaceColor',[1 1 0.5])
xlim([2 4])
set(gca,'FontSize',20)
xlabel('Distance from wall (pixels)')
ylabel('Proportion')
legend('Lights off','Lights on')
legend('boxoff')
text(3.4, 0.05,strcat({'n = '},num2str(size(allDist,1))),'FontSize',16)

shg