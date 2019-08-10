%% Kyobi Skutt-Kakaria
% 01.22.2018 - Created
% 01.23.2018 - Updated
% Harvard University
% de Bivort Lab
% this function calculates coefficient of variation for each phenotype across the data set, in this
% case its for high temperature


function out = generateLDMPlot(array)

for kk = 1:2
    
    
    for jj = 1:size(tempArray,2)
        m = cat(1,tempArray{:,jj});
        if size(m,2) == 1
            m1 = bootstrp(100,@nanmean,m);
            m2 = bootstrp(100,@(x) nanstd(x(:,1))/nanmean(x(:,1)),m);
        else
            % this treats every data point as a normal distribution with error set by earlier methods
            m1 = bootstrp(100,@(x) nanmean(normrnd(x(:,1),x(:,2))),m);
            f = @(x) nanstd(x)/nanmean(x);
            m2 = bootstrp(100,@(x) f(normrnd(x(:,1),x(:,2))),m);
        end
        datMat(kk,jj,1) = nanmean(m(:,1));
        datMat(kk,jj,2) = std(m1);
        datMat(kk,jj,3) = nanstd(m(:,1))/nanmean(m(:,1));
        datMat(kk,jj,4) = nanstd(m2);
        
        hist(m(:,1),30)
        drawnow
    end
end

cve1 = datMat(1:2,:,4)';
cv1 = datMat(1:2,:,3)';
barwitherr(cve1,cv1)
shg

out = [cv1 cve1];