function ballStruct = loadPhysiology(path,ballStruct)

if exist('imagestacks.mat') == 2
    delete('imagestacks.mat')
end

if isfield(ballStruct,'physiology')
    ballStruct = rmfield(ballStruct,'physiology');
end

% Change folder to the one containing the tif stacks

ballStruct.physiology.files = dir('*.tif');

file_name = strsplit(path,'/');


% Match up experiment names to find all the individual stacks for each
% trial

ballStruct.physiology.experimentalNum = nan(size(ballStruct.physiology.files,1),2);
for ii = 1:size(ballStruct.physiology.files,1)
    
    % Finds indexes within file name strings that contain trial number and
    % stack file number
    tifNameInd = regexp(ballStruct.physiology.files(ii).name,'\d_?\d\d\d_\d\d\d.tif\>');
    startInd = regexp(ballStruct.physiology.files(ii).name,'\d\d\d_\d\d\d.tif\>');
    endInd = regexp(ballStruct.physiology.files(ii).name,'.tif\>');
    
    % Extracts and splits trials from tif stacks 
    tempName = ballStruct.physiology.files(ii).name(startInd:endInd-1);
    ballStruct.physiology.experimentalNum(ii,:) = str2double(strsplit(tempName,'_'));
    
    % Extract file name for each tif file to be sure they are matched
    ballStruct.physiology.tifName{ii} = ballStruct.physiology.files(ii).name(1:(tifNameInd));
    
end

% Look for mismatches and remove from ballStruct
includedTifs = strcmp(ballStruct.physiology.tifName,file_name(6));
ballStruct.physiology.tifName = ballStruct.physiology.tifName(includedTifs);
ballStruct.physiology.experimentalNum = ...
    ballStruct.physiology.experimentalNum(includedTifs,:);

% obtains indexes for every trial
[ballStruct.physiology.uniqueTrials,~,ballStruct.physiology.experimentInd] = ...
    unique(ballStruct.physiology.experimentalNum(:,1));

% Define the function to fit light transition time
offsetFunc = fittype(@(a,b,k,x0,x) a + ( (b - a) ./ (1 + exp(-k.*(x - x0))) ) );


%% Initialize the master time points for all further analysis

ballStruct.physiology.masterImageTime = (-50:0.1:110);

% Assign time vector
try
    numFrames = ballStruct.internal.state{2}{1}.numberOfFrames;
catch
    try
    numFrames = ballStruct.internal.state.acq.numberOfFrames;
    catch
        numFrames = ballStruct.settings.experimentLength*ballStruct.internal.frameRate;
    end
end
frameRate = numFrames/ballStruct.settings.experimentLength;
timeVector = (1:numFrames)*(1/frameRate);
        

%% Load the tif stacks from hard disk

% make datastore to hold tif files for tall arrays
files = ballStruct.physiology.files;
experimentInd = ballStruct.physiology.experimentInd;
index = ballStruct.physiology.uniqueTrials';
index = index(index<=30);
% index = 1:15;
fullRStack = cell(0);
fullGStack = cell(0);
normalizedG = cell(0);
timeOffset = 0;

timeVect = ballStruct.physiology.masterImageTime;
bright_index =  timeVector > 60 & timeVector < 120;

ref = false;

for nn = index
    nn
    
    try
        experiment = nn;
        subFiles = files(experimentInd==experiment);
        Rstack = cell(0);
        Gstack = Rstack;
        
        for ii = 1:size(subFiles,1)
            tifStackSize = (length(imfinfo(subFiles(ii).name))/2);
            for kk = 1:tifStackSize
                subInd = kk;
                Rstack{ii}(:,:,kk) = imread(subFiles(ii).name,'Index',(subInd-1)*2+1);
                Gstack{ii}(:,:,kk) = imread(subFiles(ii).name,'Index',(subInd-1)*2+2);
                
            end
           
        end
              
        redstack = cat(3,Rstack{:});
        greenstack = cat(3,Gstack{:});
                 
       
        % Spatial smoothing
        temp1 = redstack;
        temp2 = greenstack;
        
        rMax = prctile(temp1(:),99.9);
        gMax = prctile(temp2(:),99.9);
        
        temp1(temp1>rMax) = nanmedian(temp1(:));
        temp2(temp2>gMax) = nanmedian(temp2(:));
                  
        pixelSmooth = size(temp1,2)*0.02;
        for kk = 1:min([size(redstack,3),size(greenstack,3)])
            temp1(:,:,kk) = imgaussfilt(temp1(:,:,kk),pixelSmooth);
            temp2(:,:,kk) = imgaussfilt(temp2(:,:,kk),pixelSmooth);
        end
        
        % Correct timing by internal change in background
        analysisRange = round([50 70]*frameRate);
        bgVectG =  mean(reshape(temp2,numel(temp2(:,:,1)),size(temp2,3)),1);
        vals = bgVectG(analysisRange(1):analysisRange(2));
        logFit = fit((analysisRange(1):analysisRange(2))',vals',...
            offsetFunc,...
            'Start',[min(vals), max(vals), 0.001,...
            mean(analysisRange)]);
        timeOffset(nn) = logFit.x0*(1/frameRate);
        timeVectorCorrected{nn} = timeVector - timeOffset(nn);
              
        temp1_ = double(reshape(temp1,numel(temp1(:,:,1)),size(temp1,3)));
        temp2_ = double(reshape(temp2,numel(temp2(:,:,1)),size(temp2,3)));
        
        temp1_mean = median(round(temp1_,1));
        temp2_mean = median(round(temp2_,1));
    
        temp1__ = double(temp1_) - temp1_mean;
        temp2__ = double(temp2_) - temp2_mean;
                               
        % Temporal smoothing
        lowpass_frequency = 1 * ballStruct.internal.frameRate; %hertz
        temp1__ = movmean(temp1__,lowpass_frequency,2);
        temp2__ = movmean(temp2__,lowpass_frequency,2);
        
        % Downsample to sampling rate            
        temp1__ = interp1(timeVectorCorrected{nn},temp1__',...
            timeVect)';
        temp2__ = interp1(timeVectorCorrected{nn},temp2__',...
            timeVect)';
        
        % Convert back to unsigned int
        temp1 = reshape(uint16(temp1__),size(redstack,1),size(redstack,2),size(temp1__,2));
        temp2 = reshape(uint16(temp2__),size(redstack,1),size(redstack,2),size(temp1__,2));
        
        
        % Spatial Registration
        redstack = temp1;
        greenstack = temp2;
        
        clear temp1 temp2
               
        
        
        if ~ref
            bright_index_downsample = logical(interp1(timeVectorCorrected{nn},...
                double(bright_index),timeVect));
            reference = nanmean(redstack(:,:,bright_index_downsample),3);
            ref = true;
            ballStruct.physiology.reference = reference;
        end

        subplot(1,3,2)
        hImage2 = imagesc(redstack(:,:,1));
        c = caxis;
        caxis(c);
        colormap(jet)
        subplot(1,3,3)
        hImage3 = imagesc(greenstack(:,:,1));
        c = caxis;
        caxis(c);
        colormap(jet)
        stackR = zeros(size(redstack));
        stackG = stackR;
        for ii = 1:size(redstack,3)
            try
                imgR = redstack(:,:,ii);
                imgG = greenstack(:,:,ii);
                
                usfac = 100;
                [output, Greg] = dftregistration(fft2(reference),fft2(imgR),usfac);
                stackR(:,:,ii) = abs(ifft2(Greg));
                
                
                [nr,nc]=size(imgG);
                Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
                [Nc,Nr] = meshgrid(Nc,Nr);
                Greg = fft2(imgG).*exp(1i*2*pi*(-output(3)*Nr/nr-output(4)*Nc/nc));
                Greg = Greg*exp(1i*output(2));
                stackG(:,:,ii) = abs(ifft2(Greg));
                
                
                if mod(ii,25) == 0
                    hImage2.CData = stackR(:,:,ii);
                    hImage3.CData = stackG(:,:,ii);
                    drawnow
                end
                
                
            end
            
        end
        

        redstack = stackR;
        greenstack = stackG;
     
        % Reshape to 2d
        reshapedR = double(reshape(redstack,numel(redstack(:,:,1)),size(redstack,3)));
        reshapedG = double(reshape(greenstack,numel(greenstack(:,:,1)),size(greenstack,3)));
        
        
        
        if nn == 1
            reference_green = mean(greenstack,3);
            reference_temp = reference_green(reference_green > 0);
            gmDist = fitgmdist(reference_temp,2);
            idxHi = reference_green >= max(gmDist.mu);
            ballStruct.physiology.neuronMask = reshape(idxHi,size(redstack,1),size(redstack,2));
            subplot(1,3,1)
            imagesc(ballStruct.physiology.neuronMask);
        end
         
        red = median(reshapedR(idxHi,:));
        green =  median(reshapedG(idxHi,:));
        
        if  std(red)/mean(red) > 0.25
            red = nan(1,length(red));
            green = nan(1,length(red));
        end
                
        % Normalizing function
        f = @(y1,y2) ((y1./nanmean(y1,2)) ./ (y2./nanmean(y2,2))) .*nanmean(y1,2);
        normalizedG = f(green,red);
        
        
        fullRStack = {reshape(reshapedR,[size(redstack(:,:,1)) size(reshapedR,2)])};
        fullGStack = {reshape(reshapedG,[size(greenstack(:,:,1)) size(reshapedG,2)])};
            
        meanframe = nanmean(greenstack,3);
        ballStruct.physiology.meanframe = meanframe;
        stdframe = nanstd(greenstack,[],3);
        ballStruct.physiology.stdframe = stdframe;
        
        normGStack = {reshapedG,reshapedR};
        
        if exist('imagestacks.mat') == 0
            save('imagestacks.mat','fullRStack','fullGStack','normGStack','-v7.3')
            matFileContent = matfile('imagestacks.mat','Writable',true);
        else
            matFileContent.fullRStack(1,nn) = fullRStack;
            matFileContent.fullGStack(1,nn) = fullGStack;
        end
         
        subplot(1,3,3)
        plot(red./mean(red),'red','LineWidth',1)
        hold on
        plot(green./mean(green),'green','LineWidth',1)
        plot(normalizedG./mean(normalizedG),'black','LineWidth',1)
        hold off
        pause(0.5)
        
        ballStruct.physiology.green{nn} = green;
        ballStruct.physiology.red{nn} = red;

        ballStruct.physiology.greenNormalized{nn} = normalizedG;
        ballStruct.physiology.normGStack{nn} = normGStack;


    catch
        disp(strcat({'Skipped '},num2str(nn)))
    end
end

ballStruct.physiology.timeVectorCorrected = timeVectorCorrected;
ballStruct.physiology.timeOffset = timeOffset;
ballStruct.physiology.timeVector = timeVector;

ballStruct.internal.frameRate = frameRate;


%%