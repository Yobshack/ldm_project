function turnDirs = extractTurns(flyTracks,files, lightDat)

%     load(strcat('./',files));
    flyTracks.name = files(1:end-4);
    
    flyTracks = adjustTimeStamps(flyTracks);
    
    [turnDirs,flyTracks] = determine_turns_120(flyTracks);
    if size(flyTracks.labels,2) == 6 && isstring(flyTracks.labels)
        turnDirs.tLabels = [cellstr(flyTracks.labels(2:end,1:4)),...
        num2cell(flyTracks.mazeNum)'];        
    elseif size(flyTracks.labels,2) == 6
        turnDirs.tLabels = [table2cell(flyTracks.labels(:,1:4)),...
        num2cell(flyTracks.mazeNum)'];
    elseif size(flyTracks.labels,2) == 5
        turnDirs.tLabels = [flyTracks.labels(:,1:4),...
        num2cell(flyTracks.mazeNum)'];
    end
    
    %turnDirs = calculateRelativeTime(turnDirs,lightDat);
    
    txt = strsplit(char(files),'-');
    turnDirs.dateAndTime = txt(1:5);
    
        
    turnDirs.lightDat = lightDat;
    turnDirs = calculateWallDist(turnDirs,flyTracks);

    

