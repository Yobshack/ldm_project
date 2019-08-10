function out = singleTrayAnalysis(matFile,lightDat)

load(matFile)
lights = importdata(lightDat);
turnDirs = extractTurns(flyTracks,matFile,lights);

[allDataCell, cellColNames] = parseProcessedFiles(turnDirs);

allDataCell = allDataCell(~cellfun(@isempty,allDataCell(:,1)),:);

f = @(x,y) abs(nanmean(x{y,1}(x{y,11})) - nanmean(x{y,1}(~x{y,11})));

ldm = nan(size(allDataCell,1),3);
for ii = 1:size(allDataCell,1)
    ldm(ii,:) = [f(allDataCell,ii), size(allDataCell{ii,1},1), allDataCell{ii,8}];
end

sorted = sortrows(ldm,'descend');

out = sorted;

{'LDM','Num Turns','Maze Num'}



