%% Kyobi Skutt-Kakaria
% 01.20.2018
% Harvard University
% de Bivort Lab

% This function calculates the average distance from the wall for each animal relative to the turn
% they are about to make. Smaller numbers mean closer passage to wall.

function out = wallAssym(walls,turns)

% wall prox index
f = @(x,y) nanmean([nanmean(x(y)),nanmean(x(~y))]);
wallAssymIdx = bootstrp(1000,f,walls,turns);


out = [nanmean(wallAssymIdx), nanstd(wallAssymIdx)];