threshIdx = cell2mat(allDataCell(:,13)) >= 100;

behaveMatLightLow = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,14:6:43)];
behaveMatDarkLow = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,15:6:43)];
behaveMatLightHigh = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,17:6:43)];
behaveMatDarkHigh = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,18:6:43)];

% %% Number of turns plot (whole data set)
% 
% clf
% [~,~,c] = unique(behaveMatLightLow(:,1))
% 
% histogram(cell2mat(behaveMatLightLow(:,3)),'Normalization','probability')
% hold on 
% histogram(cell2mat(behaveMatDarkLow(:,3)),'Normalization','probability')

 
 
 %% Make a distance correlation for each genotype and temperature
 
 clf
filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

effName = {'Shi','Trp','Iso'};

metricIdx = [4];

dataMatLow = nan(47,3);
dataMatHigh = nan(47,3);


 for ii = 1:47
     for jj = 1:3
         genoI = contains(behaveMatLightLow(:,1),importedData.textdata(1,1+ii),'IgnoreCase',true);
         effI = contains(behaveMatLightLow(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         if sum(Idx) < 100
             continue
         end
         dataMatLow(ii,jj) = distcorr(cell2mat(behaveMatLightLow(Idx,metricIdx)),...
             cell2mat(behaveMatDarkLow(Idx,metricIdx)));
         dataMatHigh(ii,jj) = distcorr(cell2mat(behaveMatLightHigh(Idx,metricIdx)),...
             cell2mat(behaveMatDarkHigh(Idx,metricIdx)));
     end
 end
 
 
 effSize = [(dataMatLow(:,1) - dataMatLow(:,3))./dataMatLow(:,3),...
     (dataMatHigh(:,1) - dataMatHigh(:,3))./dataMatHigh(:,3)]
 
 scatter(effSize(:,2),effSize(:,1))
 text(effSize(:,2), effSize(:,1), importedData.textdata(1,2:48),'FontSize',14,'Color','red');
 
  effSize2 = [(dataMatLow(:,2) - dataMatLow(:,3))./dataMatLow(:,3),...
     (dataMatHigh(:,2) - dataMatHigh(:,3))./dataMatHigh(:,3)]
 hold on
 scatter(effSize2(:,2),effSize2(:,1))
 text(effSize2(:,2), effSize2(:,1), importedData.textdata(1,2:48),'FontSize',14,'Color','blue');
 
xlim([-1 1])
ylim([-1 1])
%%

threshIdx = cell2mat(allDataCell(:,104)) >= 200;

behaveMatLightHigh1 = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,101:12:152)];
behaveMatLightHigh2 = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,102:12:152)];
behaveMatDarkHigh1 = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,103:12:152)];
behaveMatDarkHigh2 = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,104:12:152)];


 %% Make a distance correlation for each genotype and temperature
 
clf
filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

effName = {'Shi','Trp','Iso'};

metricIdx = [4];

dataMatLtD = nan(47,3);
dataMatLtL = nan(47,3);
dataMatDtD = nan(47,3);
sampleSize = dataMatLtD;

 for ii = 1:47
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatLightHigh1(:,1),importedData.textdata(1,1+ii),'IgnoreCase',true);
         effI = contains(behaveMatLightHigh1(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         if sum(Idx) < 60
             continue
         end
         dataMatLtD(ii,jj) = distcorr(cell2mat(behaveMatLightHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatDarkHigh2(Idx,metricIdx)));
         dataMatLtL(ii,jj) = distcorr(cell2mat(behaveMatLightHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatLightHigh2(Idx,metricIdx)));
         dataMatDtD(ii,jj) = distcorr(cell2mat(behaveMatDarkHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatDarkHigh2(Idx,metricIdx)));
         sampleSize(ii,jj) = sum(Idx);
     end
 end
 
 
 
 effSize = [(dataMatLtD(:,1) - dataMatLtD(:,3))./dataMatLtD(:,3),...
     (dataMatLtL(:,1) - dataMatLtL(:,3))./dataMatLtL(:,3),...
     (dataMatDtD(:,1) - dataMatDtD(:,3))./dataMatDtD(:,3)]
 
 scatter(effSize(:,1),effSize(:,2))
 text(effSize(:,1), effSize(:,2), importedData.textdata(1,2:48),'FontSize',4,'Color','red');
 
 effSize2 = [(dataMatLtD(:,2) - dataMatLtD(:,3))./dataMatLtD(:,3),...
     (dataMatLtL(:,2) - dataMatLtL(:,3))./dataMatLtL(:,3),...
     (dataMatDtD(:,2) - dataMatDtD(:,3))./dataMatDtD(:,3)]
 
%  hold on
%  scatter(effSize2(:,1),effSize2(:,2))
%  text(effSize2(:,1), effSize2(:,2), importedData.textdata(1,2:48),'FontSize',4,'Color','blue');
%  
% xlim([-0.8 0.8])
% ylim([-0.4 0.4])

dataMatLtDL = 2*(dataMatLtD./(dataMatLtL+dataMatDtD))
dataMat2 = (dataMatLtDL(:,1:2) - dataMatLtDL(:,3))./dataMatLtDL(:,3);

for ii = 1:length(importedData.textdata(2:end,1))
    for jj = 1:2
    idx = logical(importedData.data(ii,1:length(dataMat2))');
    dataMat3(ii,jj) = nanmean(dataMat2(idx,jj));
    end
end

bar(dataMat3)
set(gca,'XTickLabel',importedData.textdata(2:end,1),'XTickLabelRotation',45,'XTick',1:16)
shg
 %% Make a distance correlation for each genotype and temperature
 
clf
filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

effName = {'Shi','Trp','Iso'};

metricIdx = [4];

dataMatLtD = nan(47,3);
dataMatLtL = nan(47,3);
dataMatDtD = nan(47,3);
sampleSize = dataMatLtD;

 for ii = 1:47
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatLightHigh1(:,1),importedData.textdata(1,1+ii),'IgnoreCase',true);
         effI = contains(behaveMatLightHigh1(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         if sum(Idx) < 60
             continue
         end
         dataMatLtD(ii,jj) = distcorr(cell2mat(behaveMatLightHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatDarkHigh2(Idx,metricIdx)));
         dataMatLtL(ii,jj) = distcorr(cell2mat(behaveMatLightHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatLightHigh2(Idx,metricIdx)));
         dataMatDtD(ii,jj) = distcorr(cell2mat(behaveMatDarkHigh1(Idx,metricIdx)),...
             cell2mat(behaveMatDarkHigh2(Idx,metricIdx)));
         sampleSize(ii,jj) = sum(Idx);
     end
 end
 
 
 
 effSize = [(dataMatLtD(:,1) - dataMatLtD(:,3))./dataMatLtD(:,3),...
     (dataMatLtL(:,1) - dataMatLtL(:,3))./dataMatLtL(:,3),...
     (dataMatDtD(:,1) - dataMatDtD(:,3))./dataMatDtD(:,3)]
 
 scatter(effSize(:,1),effSize(:,2))
 text(effSize(:,1), effSize(:,2), importedData.textdata(1,2:48),'FontSize',4,'Color','red');
 
 effSize2 = [(dataMatLtD(:,2) - dataMatLtD(:,3))./dataMatLtD(:,3),...
     (dataMatLtL(:,2) - dataMatLtL(:,3))./dataMatLtL(:,3),...
     (dataMatDtD(:,2) - dataMatDtD(:,3))./dataMatDtD(:,3)]
 
 hold on
 scatter(effSize2(:,1),effSize2(:,2))
 text(effSize2(:,1), effSize2(:,2), importedData.textdata(1,2:48),'FontSize',4,'Color','blue');
 
xlim([-0.8 0.8])
ylim([-0.4 0.4])

%%

behaveNames = {'Turns (L)','Turns (D)','Turn Bias (L)','Turn Bias (D)',...
    'Clump (L)','Clump (D)','Switch (L)','Switch (D)','Wall Dist (L)','WallDist (D)'};
behaveNames2 = {'Turns','Turn Bias',...
    'Clump','Switch','Wall Dist'};


threshIdx = cell2mat(allDataCell(:,16)) >= 300;
behaveMatSub = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,interleave2(17:6:43,18:6:43))];


figure

blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);
nanticoke=interp1([1 52 103 164 225 256],[0 1 1; 0 .2 1; 0 0 0; 1 .1 0; 1 .9 0; 1 1 1],1:256);
nyidalur=interp1([1 128 129 256],[1 0 1; .2559 .248 .2559; .248 .2559 .248; 0 1 0],1:256);

data = behaveMatSub(contains(behaveMatSub(:,2),'Iso','IgnoreCase',true),:);
fullCorr = corr(cell2mat(data(:,3:end)),'rows','pairwise');
imagesc(fullCorr)
colormap(blanding)
caxis([-1 1])
colorbar
pbaspect([1 1 1])
title('Pooled Controls')
set(gca,'XTickLabel',behaveNames,'XTickLabelRotation',45,'YTickLabel',...
    behaveNames,'XTick',1:length(behaveNames),'YTick',1:length(behaveNames))


filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

effName = {'Shi','Trp','Iso'};

sampleSize = nan(47,3);
 corrMat = cell(0);
 corrMatPost = cell(0);
figure
compNum = 1;
 for ii = 1:47
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatSub(:,1),importedData.textdata(1,1+ii),'IgnoreCase',true);
         effI = contains(behaveMatSub(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         sampleSize(ii,jj) = sum(Idx);
         if sum(Idx) < 60
             continue
         end
         corrMat{ii,jj} = corr(cell2mat(behaveMatSub(Idx,3:end)),'rows','pairwise');


     end
     if all(sampleSize(ii,[compNum 3]) > 60)
     corrMatPost{ii} = corrMat{ii,compNum} - fullCorr;%corrMat{ii,3};
     subplot(8,6,ii)
     imagesc(corrMatPost{ii})
     caxis([-1 1])
     colormap(blanding)
     pbaspect([1 1 1])
     title(importedData.textdata(1,1+ii))
     set(gca,'XTickLabel',behaveNames2,'XTickLabelRotation',45,'YTickLabel',...
    behaveNames2,'XTick',1:2:length(behaveNames),'YTick',1:2:length(behaveNames),'FontSize',2)
     
     end
     neuron{ii} = find(importedData.data(:,ii));
 end
 
 corrMat = cell(0);
 corrMatPost = cell(0);
 clear neuron
   figure
 for ii = 1:16
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatSub(:,1),importedData.textdata(1,1+find(importedData.data(ii,1:47))),'IgnoreCase',true);
         effI = contains(behaveMatSub(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         sampleSize(ii,jj) = sum(Idx);
         if sum(Idx) < 60
             continue
         end
         corrMat{ii,jj} = corr(cell2mat(behaveMatSub(Idx,3:end)),'rows','pairwise');


     end
     if all(sampleSize(ii,[compNum 3]) > 60)
     corrMatPost{ii} = corrMat{ii,compNum} - fullCorr;%corrMat{ii,3};
     subplot(4,4,ii)
     imagesc(corrMatPost{ii})
     caxis([-1 1])
     colormap(blanding)
     pbaspect([1 1 1])
     title(importedData.textdata(ii+1,1))
     neuron{ii} = ii;
     set(gca,'XTickLabel',behaveNames2,'XTickLabelRotation',45,'YTickLabel',...
     behaveNames2,'XTick',1:2:length(behaveNames),'YTick',1:2:length(behaveNames))
     
     end
 end
 
 %%
 clf
 molaspass=interp1([1 51 102 153 204 256],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
 
 y = nan(length(corrMatPost),23);
 for ii = 1:length(corrMatPost)
     x = triu(corrMatPost{ii});
     if ~isempty(x)
        dat = x(x~=0);
        length(dat)
        dat = dat(1:2:end)
      y(ii,:) = dat;
     end
     
     ii
 end
 
 for ii = 1:length(neuron)
     n(ii) = neuron{ii}(1)
 end
 n = n(any(~isnan(y),2))
 y = y(any(~isnan(y),2),:);
 
 classes = [repmat(1,1,4) repmat(2,1,4) 1 repmat(3,1,2) 4 5 5 4 4]

shorterNames = {'P-FN1/2','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','PF-FGall','PF-FRub',...
    'PF-LCre','LPs-P',...
    'PBCap','PBInt','Ps-P1','Ps-P2'};
 
[a,b,c,d,e,f] = pca(y');
a = tsne(y,'Distance','spearman')
shg
%gscatter(z(:,2),z(:,3),n)
if size(y,1) > 16

    scatter(a(:,1),a(:,2),500,n,'MarkerFaceColor','flat')%,x(classes,:))
    text(a(:,1)-0.02*range(a(:,1)),a(:,2),shorterNames(n),'FontSize',7)
else

    scatter(a(:,2),a(:,1),500,classes,'MarkerFaceColor','flat')%,x(classes,:))
    text(a(:,2),a(:,1),shorterNames)
end
colormap(brewermap(size(y,1),'Spectral'))
xlabel('x')
ylabel('y')
zlabel('PC3')

shg

sum(sum(pdist2(y,y)))/(length(y).^2)

 %%
 blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);
 data = behaveMatSub(contains(behaveMatSub(:,2),'Shi','IgnoreCase',true),:);
fullCorr = corr(cell2mat(data(:,3:end)),'rows','pairwise');
imagesc(fullCorr)
colormap(blanding)
caxis([-1 1])
colorbar
pbaspect([1 1 1])
title('Pooled Controls')
set(gca,'XTickLabel',behaveNames,'XTickLabelRotation',45,'YTickLabel',...
    behaveNames,'XTick',1:length(behaveNames),'YTick',1:length(behaveNames))
 
x = tsne(cell2mat(data(:,3:end)));
y = data(:,1)

 gscatter(x(:,1),x(:,2),y)
 
 %%

behaveNames = {'Turns (L)','Turns (D)','Turn Bias (L)','Turn Bias (D)',...
    'Clump (L)','Clump (D)','Switch (L)','Switch (D)','Wall Dist (L)','WallDist (D)'};
behaveNames2 = {'Turns','Turn Bias',...
    'Clump','Switch','Wall Dist'};


threshIdx = cell2mat(allDataCell(:,16)) >= 500;
behaveMatSub = [allDataCell(threshIdx,3:4) allDataCell(threshIdx,interleave2(17:6:43,18:6:43))];


figure


nanticoke=interp1([1 52 103 164 225 256],[0 1 1; 0 .2 1; 0 0 0; 1 .1 0; 1 .9 0; 1 1 1],1:256);
nyidalur=interp1([1 128 129 256],[1 0 1; .2559 .248 .2559; .248 .2559 .248; 0 1 0],1:256);

data = behaveMatSub(contains(behaveMatSub(:,2),'Iso','IgnoreCase',true),:);
fullCorr = corr(cell2mat(data(:,3:end)),'rows','pairwise');
imagesc(fullCorr)
colormap(blanding)
caxis([-1 1])
colorbar
pbaspect([1 1 1])
title('Pooled Controls')
set(gca,'XTickLabel',behaveNames,'XTickLabelRotation',45,'YTickLabel',...
    behaveNames,'XTick',1:length(behaveNames),'YTick',1:length(behaveNames))


filepath = '/Users/kyobikakaria/Desktop/Data/split_screen/rubinAnnotation.csv';
importedData = importdata(filepath);

effName = {'Shi','Trp','Iso'};

sampleSize = nan(47,3);
 corrMat = cell(0);
 corrMatPost = cell(0);
 metricMat = nan(47,30);
compNum = 1;
 for ii = 1:47
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatSub(:,1),importedData.textdata(1,1+ii),'IgnoreCase',true);
         effI = contains(behaveMatSub(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         sampleSize(ii,jj) = sum(Idx);
         if sum(Idx) < 60
             continue
         end
         data = cell2mat(behaveMatSub(Idx,3:end));
         corrMat{ii,jj} = [nanmean(data,1), nanstd(data,[],1),...
             nanmean(data(:,1:2:end)-data(:,2:2:end)), nanstd(data(:,1:2:end)-data(:,2:2:end),[],1)];


     end
     if all(sampleSize(ii,[compNum 3]) > 60)
     metricMat(ii,:) = (corrMat{ii,compNum} - corrMat{ii,3} + 0.1)./(corrMat{ii,3}+0.1);

     end
     neuron{ii} = find(importedData.data(:,ii));
 end
 
 metricMat = nan(16,30);
 corrMat = cell(0);
 corrMatPost = cell(0);

 for ii = 1:16
     genotypes(ii) = importedData.textdata(1,1+ii);
     for jj = 1:3
         genoI = contains(behaveMatSub(:,1),importedData.textdata(1,1+find(importedData.data(ii,1:47))),'IgnoreCase',true);
         effI = contains(behaveMatSub(:,2),effName(jj),'IgnoreCase',true);
         Idx = genoI & effI;
         sampleSize(ii,jj) = sum(Idx);
         if sum(Idx) < 60
             continue
         end
         data = cell2mat(behaveMatSub(Idx,3:end));
         corrMat{ii,jj} = [nanmean(data,1), nanstd(data,[],1),...
             nanmean(data(:,1:2:end)-data(:,2:2:end)), nanstd(data(:,1:2:end)-data(:,2:2:end),[],1)];


     end
     if all(sampleSize(ii,[compNum 3]) > 60)
     metricMat(ii,:) = (corrMat{ii,compNum} - corrMat{ii,3})./corrMat{ii,3};

     end
 end
 
 %%
 clf
 molaspass=interp1([1 51 102 153 204 256],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
 

 for ii = 1:length(neuron)
     n(ii) = neuron{ii}(1)
 end
 n = n(any(~isnan(y),2))
  y = metricMat
 y = y(any(~isnan(y),2),:);
 

 classes = [repmat(1,1,4) repmat(2,1,4) 1 repmat(3,1,2) 4 5 5 4 4]

shorterNames = {'P-FN1/2','P-FN3','P-FN4','P-FN5','P-EN','P-EG','E-PG','E-PG9','PF-FGall','PF-FRub',...
    'PF-LCre','LPs-P',...
    'PBCap','PBInt','Ps-P1','Ps-P2'};
 
[a,b,c,d,e,f] = pca(y');
%a = tsne(y)
shg
%gscatter(z(:,2),z(:,3),n)
if size(y,1) > 16

    scatter(a(:,1),a(:,2),500,n,'MarkerFaceColor','flat')%,x(classes,:))
    text(a(:,1)-0.02*range(a(:,1)),a(:,2),shorterNames(n),'FontSize',7)
else

    scatter(a(:,2),a(:,1),500,classes,'MarkerFaceColor','flat')%,x(classes,:))
    text(a(:,2),a(:,1),shorterNames)
end
colormap(brewermap(size(y,1),'Spectral'))
xlabel('x')
ylabel('y')
zlabel('PC3')

shg

sum(sum(pdist2(y,y)))/(length(y).^2)


