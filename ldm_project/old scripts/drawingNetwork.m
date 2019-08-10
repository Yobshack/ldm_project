function drawingNetwork(input)


filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/split_Rubin_AnnotationCons/Sheet 2-Table 1.csv';
importedData = importfile1(filepath);
array = table2array(importedData(2:end,2:end));

names = importedData(:,1);

shorterNeuropilNames = {'PB','FB','EB','NO','Rubus','Crepine','Gall','LAL','PS'};

shorterNames = {'P-FN1/2','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','PF-FGall','PF-FRub',...
    'PF-LCre','LPs-P',...
    'PBCap','PBInt','Ps-P1','Ps-P2','EB','LAL','Rub','Cre','PS','Gall'};

PBOutput = [1 2 3 4 5 6 9 10 11];
PBInput = [7 8 12 15 16];
PBIntrinsic = [13 14];
projectionNeuropil = [17:22];

reorder = [PBInput, PBIntrinsic, PBOutput];

array = array(reorder,reorder);

G = digraph(array~=0);

p = plot(G,'Layout','circle');
p.ArrowSize = 5;


kotor=interp1([1 84 128 129 173 256],[0 0 1; 0 1 1; .248 .2559 .2559; .248 .248 .2559; 1 1 0; 1 0 0],1:256);
colormap(kotor)



p.NodeCData = input(reorder);

p.EdgeColor = [0 0 0];
p.LineWidth = 0.5;
p.EdgeAlpha = 0.8;
% 
% xlim([-1 1.5])
% ylim([-2 1.5])

p.NodeLabel = shorterNames(reorder);
p.MarkerSize = 5;
set(gca,'FontSize',8)

maxRange = max([max(input) abs(min(input))]);

caxis([-maxRange maxRange])

colorbar

prob = normcdf(input(reorder)./inputE(reorder));
prob(prob>0.5) = 1-prob(prob>0.5);

pbaspect([1 1 1])

for ii = 1:length(prob)
    if prob(ii) <= 0.05 && prob(ii) > 0.01
        p.NodeLabel(ii) = strcat(p.NodeLabel(ii),'*');
    elseif prob(ii) <= 0.01 && prob(ii) > 0.001
        p.NodeLabel(ii) = strcat(p.NodeLabel(ii),'**');
    elseif prob(ii) <= 0.001
        p.NodeLabel(ii) = strcat(p.NodeLabel(ii),'***');
    end
end


shg