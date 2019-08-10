function out = aggregateStatistics(inputStruct,indexes);

subStruct = inputStruct(indexes);
turnBiasesLowLight = extractTB(subStruct);

subStruct = inputStruct.lowTemp.dark(indexes);
turnBiasesLowDark = extractTB(subStruct);

subStruct = inputStruct.highTemp.light(indexes);
turnBiasesHighLight = extractTB(subStruct);

subStruct = inputStruct.highTemp.dark(indexes);
turnBiasesHighDark = extractTB(subStruct);

tbExp = [turnBiasesLowLight turnBiasesLowDark turnBiasesHighLight turnBiasesHighDark];

% extracts mads in same order as above
mads = calculateMAD(tbExp);
ldms = calculateLDM(tbExp);



