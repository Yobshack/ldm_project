function ballStruct = analyzeBehavior(ballStruct,path)

%% Analyze behavior and plot activity for individual


% initialize ball motion variables
ballStruct.behavior.smoothedYaw = cell(0);
ballStruct.behavior.smoothedPitch = ballStruct.behavior.smoothedYaw;
ballStruct.behavior.smoothedRoll = ballStruct.behavior.smoothedYaw;
for nn = 1:size(ballStruct.behavior.unionExperiments,1)
    % Obtain frame rate
%     if strcmp(ballStruct.settings.acqMode,'Plane')
        frameRate = ballStruct.internal.frameRate;
%     else
%         frameRate = nanmean([ballStruct.internal.frameRateVol{nn}]);
%     end
    nn
    try
        % load .dat file for an individual trial
        videoIdx = ballStruct.behavior.videoIndex(nn);
        fileName = fullfile(path,'videos',...
            ballStruct.behavior.files(videoIdx).name);
        ballStruct.behavior.datFile{nn} = importdata(fileName);
        
        % subtract estimates of the offset between behavior and physiology data
        ballStruct.behavior.behaviorTime{nn} = ballStruct.behavior.datFile{nn}(:,22)./1000 -...
            ballStruct.internal.scanImageDelay(ballStruct.behavior.imageIndex(nn)) - ballStruct.physiology.timeOffset(nn);
        
        % determine smoothing for behavior traces
        smoothTime = 3; % time in seconds to smooth, not sure if this is a good amount of time? Reasonable Q for Ben.
        ballFPS = 1/nanmean(diff(ballStruct.behavior.behaviorTime{nn}));
        ballStruct.behavior.smoothNum(nn) = smoothTime*ballFPS/frameRate;
        
        dat1 = ballStruct.behavior.datFile{nn}(:,4);
        %         dat1(abs(dat1) > nanstd(dat1)*3) = nan;
        %         dat1 = fillmissing(dat1,'linear');
        %         ballStruct.behavior.smoothedRoll{nn} = smoothdata(dat1,1,'movmean',ballStruct.behavior.smoothNum(nn));
        ballStruct.behavior.smoothedRoll{nn} = movmean(dat1,ballFPS);
        
        dat2 = ballStruct.behavior.datFile{nn}(:,3);
        %         dat2(abs(dat2) > nanstd(dat2)*3) = nan;
%         dat2 = fillmissing(dat2,'linear');
        ballStruct.behavior.smoothedYaw{nn} = movmean(dat2,ballFPS);
        
        dat3 = ballStruct.behavior.datFile{nn}(:,2);
%         dat3(abs(dat3) > nanstd(dat3)*3) = nan;
%         dat3 = fillmissing(dat3,'linear');
        ballStruct.behavior.smoothedPitch{nn} = movmean(dat3,ballFPS);
        
    catch
        disp('Failed at time standardization and matching')
        
    end

end