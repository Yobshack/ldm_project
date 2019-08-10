%% this function is used to calculate the number of turns made by an animal

function out = turnrateTTA(turns,timeBin)

turnrate = length(turns)/timeBin

out = [turnrate];
