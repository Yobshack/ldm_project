%% Kyobi Skutt-Kakaria
% 01.20.2018
% Harvard University
% de Bivort Lab

% This function calculates the average distance from the wall for each animal relative to the turn
% they are about to make. Smaller numbers mean closer passage to wall.

function out = wallDist(walls,turns)

% % % flip left turns to right turns
walls(~turns) = 1-walls(~turns);

% extract abolute wall proximity rather than 

% return average distance
out1 = nanmean(walls);

% calculate error by bootstrapping
out2 = std(bootstrp(100,@nanmean,walls));

out = [out1 out2];
