function LDM_calcium_imaging(path,ballStruct)

matFileName = ballStruct.matFileName;

ballStruct.combined.masterImageTime = (-50:0.1:110);

% Physiology Analysis
cd(strcat(path,'/tifs/'))
ballStruct = loadPhysiology(path,ballStruct);
commandLineProgressBar(1,5,'$')

% Behavior analysis 
ballStruct = loadBehavior(path,ballStruct);
commandLineProgressBar(3,5,'$')

% Analyze behavioral data
ballStruct = analyzeBehavior(ballStruct,path);
commandLineProgressBar(4,5,'$')

% Combine behavior and physiology
ballStruct = combineBehaviorPhys(ballStruct);

save(matFileName,'ballStruct','-v7.3')
