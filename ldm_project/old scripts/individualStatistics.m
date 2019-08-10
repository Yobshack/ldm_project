function out = individualStatistics(varargin)

turns = varargin{1};
normTime = varargin{2};
time = varargin{3};
if length(varargin) == 4
walls = varargin{4};
end

out.TB = nanmean(turns);
out.numTurns = length(turns);
n = out.numTurns;
p = out.TB;
out.bse = sqrt(n*p*(1-p))/n;


% calculate a transition kernel for each individual
[kernel,subTurns,bse2] = calculateKernels(normTime,turns);

out.kernel = kernel';
out.kernelTurns = subTurns';
out.bse2 = bse2;


% calculate mutual information for each animal (or auto correlation?)
[a,b,c] = autocorr(turns,10);
out.autoCorrelation = a;
out.autoCorrelationCI = c;

% clumpiness calculation
out.clump = mad(diff(time),1)./(range(time)/60);

% switchiness calculation
out.switch = nanmean((turns(1:(end-1)) & turns(2:end)))/nanmean(turns);

% intra-individual variation

firstHalf = time(1) + range(time)/2;
secondHalf = time(end) - range(time)/2;

turns1st = turns(time < firstHalf);
turns2nd = turns(time >= secondHalf);

tb1st = nanmean(turns1st);
tb2nd = nanmean(turns2nd);

out.tbSplit = [tb1st tb2nd];



%calculate predictive index of wall following

if exists(walls)
logicTable = [walls(~isnan(walls)) turns(~isnan(walls))];
out.wallPred = sum(logicTable(:,1) == logicTable(:,2))/size(logicTable,1);
end





