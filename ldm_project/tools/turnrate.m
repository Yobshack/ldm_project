%% this function is used to calculate the number of turns made by an animal

function out = turnrateTTA(turns,tStamps,lightDat)

% I will make a null distribution in which i randomly select time stamps and sum them until i exceed
% the total time of the experiment by bootstrapping which should give me a standard error on the
% estimate

diffs = [0; diff(tStamps)];


% generate transition time vector
tt = [0; cumsum(lightDat(1:(end-1),1))];

% subset by transitions that change light condition
tt0 = find([0; diff(lightDat(:,2))] ~= 0);
tt = [0; tt(tt0)];

% return turns that happen on either side of a time transition
ti = zeros(length(tStamps)-1,2);
for ii = 1:(length(tStamps)-1)
    ti(ii,:) = [sum(tStamps(ii) > tt),sum(tStamps(ii+1) > tt)];
end   

idx = ti(:,1) == ti(:,2);
diffs = diffs(idx);
turns = turns(idx);

for ii = 1:100
tTime = 0;
turnV = 1;
count = 1;
while tTime < sum(lightDat(unique(ti),1))
    idx = randi(length(turns),1);
    turnV(count) = turns(idx);
    tTime = tTime + diffs(idx);
    count = count+1;
end

reSamp(ii) = length(turnV)/sum(lightDat(unique(ti),1));

end

out = [nanmean(reSamp) nanstd(reSamp)];
