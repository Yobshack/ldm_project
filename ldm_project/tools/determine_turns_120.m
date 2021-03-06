function [out,out2] = determine_turns_120(flyTracks)
turns = flyTracks.rightTurns';
tStamps = flyTracks.tStamps;
nFlies = flyTracks.nFlies;
mazeOri = flyTracks.mazeOri;
if isfield(flyTracks,'lightBlocks')
    lightBlocks = flyTracks.lightBlocks;
end
turnIdx = ~isnan(turns);

if size(turnIdx,2) > size(turnIdx,1)
    turnIdx = turnIdx';
    turns = turns';
end

for i=1:nFlies
    idx = turnIdx(:,i);
    tSeq=turns(idx,i);
    timeSeq = turnIdx;
    tSeq=diff(tSeq);
    tt1 = tStamps(idx); tTime{i} = tt1(2:end);
    if isfield(flyTracks,'lightBlocks')
        tt2 = lightBlocks(idx,:); tLights{i} = tt2(2:end,:);
    end
    if mazeOri(i)
        tSequence{i}=tSeq==1|tSeq==-2;
    elseif ~mazeOri(i)
        tSequence{i}=tSeq==-1|tSeq==2;
    end
    
    tIdx{i} = find(idx(2:end));

end



out.tSequence = tSequence';

out.tTime = tTime';
if isfield(flyTracks,'lightBlocks')
    out.tLights = tLights';
end

out.tIdx = tIdx';

flyTracks.mazeNum = 1:120;

flyTracks.turnIdx = turnIdx;
out2 = flyTracks;
