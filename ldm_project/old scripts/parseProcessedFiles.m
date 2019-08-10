%% Kyobi Skutt-Kakaria
% 01.21.2018 - created
% 01.22.2018 - updated
% Harvard University
% de Bivort Lab
% this function combines all the processed mat files to create a single cell array with raw data from
% every fly


function [out, out2] = parseProcessedFiles(files)

tmpC = cell(0);
for ii = 1:length(files)
    if iscell(files)
        load(files{ii})
    elseif isstruct(files)
        load(files(ii).name)
    end
    ii
    % need to heal a few of the labels files because the number of roi's was not 120 on some trays
    if size(turnDirs.tLabels,1) > size(turnDirs.tSequence,1)
        turnDirs.tLabels = turnDirs.tLabels(1:size(turnDirs.tSequence,1),:);
    elseif size(turnDirs.tLabels,1) < size(turnDirs.tSequence,1)
        turnDirs.tLabels = [turnDirs.tLabels;repmat(turnDirs.tLabels(end,:),...
            size(turnDirs.tSequence,1)-size(turnDirs.tLabels,1))];
    end     
    
    if istable(turnDirs.tLabels)
        turnDirs.tLabels = table2cell(turnDirs.tLabels);
    end
        
    
    
    dataCell = cell(size(turnDirs.tSequence,1),12);
    for kk = 1:size(turnDirs.tSequence,1)
        if length(turnDirs.avgWallDist) < kk
            continue
        end
        % separate geno and effector
        labs = turnDirs.tLabels(kk,:);
        if contains(labs(1),'_')
        labs(6) = cell(1);
        labsSplit = strsplit(char(labs(1)),'_');
        labs = [labsSplit(1) upper(labsSplit(2)) labs(2:(end-1))];
        else
            labs = [labs(1) {''} labs(2:end)];
        end
        dataCell(kk,:) = [turnDirs.tSequence(kk),turnDirs.tTime(kk),labs,...
            {turnDirs.dateAndTime},turnDirs.avgWallDist(kk),...
            {logical(turnDirs.avgLight{kk})},turnDirs.proArea(kk)];
    end
    tmpC{ii} = dataCell;
end

cellColNames = {'turnSequence','turnTimes','genotype','effector','sex',...
    'treatment','box','maze','dateAndTime','wallDistance','lightStatus',...
    'totalMazeArea','relativeTime','neuronName','shorterNeuronName'};


% concatenate all data into the same cell array
allDataCell = cat(1,tmpC{:});

% take only non-empty cells
allDataCell = allDataCell(~sum(cellfun(@isempty,allDataCell(:,1:2)),2) > 0,:);

% filter for proper maze size
a = cell2mat(allDataCell(:,12));
idx = a > (median(a) - 2*mad(a)) & a < (median(a) + 2*mad(a));


cellArray = allDataCell(idx,:);


out = cellArray;
out2 = cellColNames;