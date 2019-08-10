function out = extractTB(subStruct)

turnBiases = nan(length(subStruct),1);

for ii = 1:length(turnBiases)
    if ~isempty(subStruct{ii})
        turnBiases(ii) = subStruct{ii}.TB;
    else
        turnBiases(ii) = nan;
    end
end

out = turnBiases;