function out = calculateLDM(turnBiases)


% Calculate ldm of turn bias
ldmLow = turnBiases(:,2) - turnBiases(:,1);
ldmHigh = turnBiases(:,4) - turnBiases(:,3);

reMad = nan(100,size(madTB,2));
for ii = 1:100
    row = randi(size(turnBiases,1),size(turnBiases,1),1);
    reMad(ii,:) = mad(turnBiases(row,:),1);
end

out.madTB = madTB;
out.madSEM =  nanstd(reMad);