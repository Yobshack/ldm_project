function [out1,out2,out3]=calculateKernels(time,turn)

% get kernels for each transition for each fly
turns = turn;
[a,b] = discretize(time,[0 1 4 9 21 60]);

tb = nan(length(b)-1,1);
nTurns = tb;
for jj = 1:(length(b)-1)
    tb(jj) = nanmean(turns(a==jj));
    nTurns(jj) = length(turns(a==jj));
    bse(jj) = sqrt(nTurns(jj)*tb(jj)*(1-tb(jj)))/nTurns(jj);
end

out1 = tb;
out2 = nTurns;
out3 = bse;