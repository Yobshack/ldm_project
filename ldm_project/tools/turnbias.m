%% this function calculates the turn bias of individual animals

function out = turnbias(turns)

% calculate right turn probability
out1 = nanmean(turns);

% calculate binomial sampling error
out2 = sqrt((out1*(1-out1))/length(turns));

out = [out1 out2];