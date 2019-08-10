%% Kyobi Skutt-Kakaria
% 01.20.2018
% Harvard University
% de Bivort Lab
% this function calculates the clumpiness score for an individual
% This function calculates the distribution of observed inter turn intervals compared to what would
% be expected given a metronomic turn rate, as in turns at a regular interval


function  out2 = clumpiness(tStamps,lightDat)

% set a turn filter to keep the noise away
turnfilter = 100;

% generate transition time vector
tt = [0; cumsum(lightDat(1:(end-1),1))];

% subset by transitions that change light condition
tt0 = find([1; diff(lightDat(:,2))] ~= 0);
tt = tt(tt0);

% return turns that happen on either side of a time transition
ti = zeros(length(tStamps)-1,2);
for ii = 1:(length(tStamps)-1)
    ti(ii,:) = [sum(tStamps(ii) > tt),sum(tStamps(ii+1) > tt)];
end
    
% take differences of tStamps
diffs = diff(tStamps);

% remove transitions tStamps;
diffs(~(ti(:,1) == ti(:,2))) = [];

% calculate MAD of inter-turn-interval - will try to find a way to display both the median and mean
% MAD maybe
if length(diffs) > turnfilter;
%m1 = mad(diffs,0); % mean MAD
%m1b = bootstrp(100,@(x) mad(x,0),diffs);
m2 = mad(diffs,1); % median MAD
m2b = bootstrp(100,@(x) mad(x,1),diffs);

% normalization factor - returned as turns per second
time = sum(lightDat(tt0(unique(ti)),1));
n1 = time/(length(diffs)+1);

% clumpiness score - mad of turn time normalized by the number of turns in that time window, should
% make the stat invariant with respect to num turns.

%out1 = [m1/n1 std(m1b/n1)];
out2 = [m2/n1 std(m2b/n1)];
else
    out2 = [nan nan];
end

