function ballStruct = combineBehaviorPhys(ballStruct)

if isfield(ballStruct,'combined')
    ballStruct = rmfield(ballStruct,'combined');
end


%% Sample behavior and physiology to the same rate, set by the limits of the physiology

timeVect = ballStruct.physiology.masterImageTime;

% Loop through trials to build a matrix of GCaMP values

ballStruct.combined.GCaMP = nan(length(ballStruct.behavior.smoothedPitch),...
    length(timeVect));
ballStruct.combined.yaw = ballStruct.combined.GCaMP;
ballStruct.combined.pitch = ballStruct.combined.GCaMP;
ballStruct.combined.roll = ballStruct.combined.GCaMP;


% physiology data
ballStruct.combined.GCaMP = cell2mat(ballStruct.physiology.greenNormalized');

for nn = 1:size(ballStruct.combined.GCaMP,1)

       
        try

            
            % behavior data
            yawData = cumsum(ballStruct.behavior.smoothedYaw{nn});
            pitchData = cumsum(ballStruct.behavior.smoothedPitch{nn});
            rollData = cumsum(ballStruct.behavior.smoothedRoll{nn});
            behaviorTimePts = ballStruct.behavior.behaviorTime{nn};
            
            % Sample the GCamp and behavior at the same rate
 
            
            ballStruct.combined.yaw(ballStruct.behavior.imageIndex(nn),:) = [0 diff(interp1(behaviorTimePts,yawData,...
                timeVect))];
            ballStruct.combined.pitch(ballStruct.behavior.imageIndex(nn),:) = [0 diff(interp1(behaviorTimePts,pitchData,...
                timeVect))];
            ballStruct.combined.roll(ballStruct.behavior.imageIndex(nn),:) = [0 diff(interp1(behaviorTimePts,rollData,...
                timeVect))];
            
           stackG = ballStruct.physiology.normGStack{nn}{1};
           stackR = ballStruct.physiology.normGStack{nn}{2};
           mask = ballStruct.physiology.neuronMask;
           dims = size(ballStruct.physiology.reference);
           stack = stackG./stackR;
           corrmap = reshape(corr(stack',ballStruct.combined.pitch(nn,:)'),dims(1),dims(2));
           ballStruct.combined.corrmap{nn} = corrmap;
           
           
           
%            imagesc(corrmap)
           pause(0.1)
            
         
            
        catch
            disp('Couldnt process behavior, has it been analyzed?')
        end
        
    
end

skalafell=interp1([1 128 129 256],[1 0 1; 1 .9921 1; .9921 1 .9921; 0 1 0]*0.8,1:256);

corrmap = mean(cat(3,ballStruct.combined.corrmap{:}),3);
corrmap(~mask) = nan;
im = imagesc(corrmap);
im.AlphaData = mask;
colormap(skalafell)
colorbar
drawnow

ballStruct.physiology = rmfield(ballStruct.physiology,'normGStack');

