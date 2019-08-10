x = cat(1,lightBias{:});
y = cat(1,darkBias{:});

xT = cat(1,lightTurns{:});
yT = cat(1,darkTurns{:});

ldm = x(:,2) - y(:,2);
goodIdx = ~isnan(ldm);

se = sqrt(xT(:,2).*x(:,2).*(1-x(:,2))./xT(:,2).^2 + xT(:,2).*x(:,2).*(1-x(:,2))./yT(:,2).^2);

outMat = [ldm se];

outMat = outMat(goodIdx,:);

nNames = cat(2,neuronName{:});
gNames = cat(2,genoName{:});
eNames = cat(2,effectorName{:});

nNames = nNames(goodIdx)';
gNames = gNames(goodIdx)';
eNames = eNames(goodIdx)';

[uniqueNNames,~,nI] = unique(nNames);
[uniqueGNames,~,gI] = unique(gNames);
[uniqueENames,~,eI] = unique(eNames);

%% First try of this - alot of these predictors are nested
predictionMat = zeros(size(outMat,1),...
    length(uniqueNNames)+length(uniqueGNames)+length(uniqueENames));
for ii = 1:length(uniqueNNames)
    predictionMat(nI==ii,ii) = 1;
end

for ii = 1:length(uniqueGNames)
    predictionMat(gI==ii,ii+length(uniqueNNames)) = 1;
end

for ii = 1:length(uniqueENames)
    predictionMat(eI==ii,ii+length(uniqueNNames)+length(uniqueGNames)) = 1;
end

%% Just for neurons with 1 or 0 in all cols

predictionMat = zeros(size(outMat,1),...
    length(uniqueNNames)+length(uniqueENames));

for ii = 1:length(uniqueNNames)
    predictionMat(nI==ii,ii) = 1;
end

for ii = 1:length(uniqueENames)
    predictionMat(eI==ii,ii+length(uniqueNNames)) = 1;
end

%% Just for gal4's with 1 or 0 in all cols

predictionMat = zeros(size(outMat,1),...
    length(uniqueGNames)+length(uniqueENames));

for ii = 1:length(uniqueGNames)
    predictionMat(gI==ii,ii) = 1;
end

for ii = 1:length(uniqueENames)
    predictionMat(eI==ii,ii+length(uniqueGNames)) = 1;
end

%% %% Just for gal4's with 0 in Iso col

predictionMat = zeros(size(outMat,1),...
    length(uniqueGNames)+length(uniqueENames));

for ii = 1:length(uniqueGNames)
    predictionMat(gI==ii,ii) = 1;
end

for ii = 1:length(uniqueENames)
    predictionMat(eI==ii,ii+length(uniqueGNames)) = 1;
end

predictionMat(:,length(uniqueGNames)+1) = [];

%%
path = "/Users/kyobikakaria/Documents/R/";
dlmwrite(strcat(path,'predictionMat.csv'),predictionMat)
dlmwrite(strcat(path,'ldmForStan.csv'),outMat)


