function ballStruct = loadBehavior(path,ballStruct)
%% Behavior analysis
if isfield(ballStruct,'behavior')
    ballStruct = rmfield(ballStruct,'behavior');
end

file_name = strsplit(path,'/');

% Set path to analyzed video dat files
ballStruct.behavior.files = dir(fullfile(path,'videos','*.dat'));
% this is to filter for only FicTrac V2
index = arrayfun(@(x) length(x.name) > 65, ballStruct.behavior.files);
ballStruct.behavior.files = ballStruct.behavior.files(index);

% This is needed to match behavior files to their proper tif stacks
if length(ballStruct.behavior.files(1).name) < 65
    ballStruct.behavior.videoNums = nan(length(ballStruct.behavior.files),1);
    for ii = 1:length(ballStruct.behavior.files)
        % find the file number for each file
        startInd = regexp(ballStruct.behavior.files(ii).name,'\d\d?.dat');
        endInd = regexp(ballStruct.behavior.files(ii).name,'.dat');
        ballStruct.behavior.videoNums(ii) = str2double(ballStruct.behavior.files(ii).name...
            (startInd:(endInd-1)));
        % extract the file name from the video file
        ballStruct.behavior.videoNames{ii} = ballStruct.behavior.files(ii).name(1:(startInd-2));
    end
    
else
    
    % This is for newer FicTrac v2
    for ii = 1:length(ballStruct.behavior.files)
        % find the file number for each file
        startInd = regexp(ballStruct.behavior.files(ii).name,'\d\d?-');
        endInd = regexp(ballStruct.behavior.files(ii).name,'-');
        ballStruct.behavior.videoNums(ii) = str2double(ballStruct.behavior.files(ii).name...
            (startInd:(endInd-1)));
        % extract the file name from the video file
        ballStruct.behavior.videoNames{ii} = ballStruct.behavior.files(ii).name(1:(startInd-2));
    end
    
end

% Look for mismatches and remove from ballStruct
includedVideos = strcmp(ballStruct.behavior.videoNames,file_name(6));
ballStruct.behavior.videoNames = ballStruct.behavior.videoNames(includedVideos);
ballStruct.behavior.videoNums = ballStruct.behavior.videoNums(includedVideos);

% find where tif and video file names match
[unionExperiments,imageIndex,videoIndex] = intersect(ballStruct.physiology.uniqueTrials,...
    ballStruct.behavior.videoNums);

ballStruct.behavior.unionExperiments = unionExperiments;
ballStruct.behavior.imageIndex = imageIndex;
ballStruct.behavior.videoIndex = videoIndex;