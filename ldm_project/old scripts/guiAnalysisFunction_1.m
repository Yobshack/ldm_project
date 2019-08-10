function out = guiAnalysisFunction_1(handles)


allGenos = handles.importedAnno.textdata(1,2:end);
% the purpose of this section is to return the indexes of flies within each grp
genoVects = cell(0);
for jj = 1:length(handles.types)
    if length(handles.types{jj}) > 7
        % this returns all the genotypes contained within a neuron ID to query the dataset with
        genos = allGenos(logical(handles.importedAnno.data(strcmp(handles.types(1),...
            handles.importedAnno.textdata(2:end,1)),:)));
    else
        genos = handles.types(jj);
    end
    if length(handles.types) > 1
        logInd = cellfun(@(x) findFlys(x,handles),genos,'UniformOutput',false);
        genoVects{jj} = any(cat(2,logInd{:}),2);
    else
        genoVects{jj} = cellfun(@(x) findFlys(x,handles),genos,'UniformOutput',false);
    end
end





switch b % Analysis grp
    case 1 % Shibire versus Control
        
        
    case 2 % Trp versus Control
    case 3 % Shibire
       
    case 4 % Trp
    case 5 % Control
end


switch c % Which turns
    case 1 % All Turns
    case 2 % Light Turns
    case 3 % Dark Turns
end

switch d % Z score option
    case 1 % Raw
    case 2 % Z score
end

switch g % temperature
    case 1 % Low
    case 2 % High
    case 3 % Low versus High
end



