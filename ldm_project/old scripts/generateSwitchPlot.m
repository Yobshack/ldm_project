%% Kyobi Skutt-Kakaria
% 02.07.2018 - Created
% 02.07.2018 - updated

% this function generates a plot of switchiness over inter-turn delay

function out = generateSwitchPlot(array,lightDat)
timeVect = 2:0.5:20; timeBin = 5;
mean1 = nan(size(array,1),length(timeVect)-1);
mean2 = mean1;
for hh = 1:size(array,1)
nums = randi(size(array,1),10,1);
% for bb = 1:10
%     hh = nums(bb);
tStamps = array{hh,2};

idx = tStamps >= 0*60 & tStamps < 80*60;

if sum(idx) < 20
    continue
end

turns = array{hh,1}(idx,:);
tStamps = tStamps(idx);

% generate transition time vector
tt = [0; cumsum(lightDat(1:(end-1),1))];

% subset by transitions that change light condition
tt0 = find([0; diff(lightDat(:,2))] ~= 0);
tt = [0; tt(tt0)];

% return turns that happen on either side of a time transition
ti = zeros(length(tStamps)-1,2);
for ii = 1:(length(tStamps)-1)
    ti(ii,:) = [sum(tStamps(ii) > tt),sum(tStamps(ii+1) > tt)];
end   

% number of right turns following right turns
rrturns = turns(1:(end-1)) + turns(2:end)==1;

% remove turns found across transitions;
rrturns(ti(:,1) ~= ti(:,2)) = [];
turns([true; ti(:,1) ~= ti(:,2)]) = [];

% diff of t stamps
distances = diff(tStamps);
distances(ti(:,1) ~= ti(:,2)) = [];

% set light and dark turn idx
ldlogical = array{hh,11}(idx);
idxL = ldlogical(2:end);
idxL(ti(:,1) ~= ti(:,2)) = [];

% normalization prob
probR = nanmean(turns(idxL));
nfL = 2.*probR.*(1-probR);
probR = nanmean(turns(~idxL));
nfD = 2.*probR.*(1-probR);

switchVect1 = nan(length(timeVect)-1,1);
switchVect2 = nan(length(timeVect)-1,1);
for ii = 1:(length(timeVect)-1)
    idx1 = any(distances(idxL)>=(timeVect(ii)-timeBin) & ...
        distances(idxL)<=(timeVect(ii+1)+timeBin),2);
    idx2 = any(distances(~idxL)>=(timeVect(ii)-timeBin) & ...
        distances(~idxL)<=(timeVect(ii+1)+timeBin),2);
    if sum(idx1) > 20
        rtt = rrturns(idxL);
        trn1 = rtt(idx1);
        rtt = rrturns(~idxL);
        trn2 = rtt(idx2);
        % number of right turns following right turns
        switchVect1(ii) = nanmean(trn1)/nfL;
        switchVect2(ii) = nanmean(trn2)/nfD;
    end
end

mean1(hh,:) = switchVect1;
mean2(hh,:) = switchVect2;

end

%%
clf
hold on

idx = all(isnan(mean1) & isnan(mean2),1);
mean1 = mean1(:,~idx);
mean2 = mean2(:,~idx);

xVect = timeVect(1:(end-1))+(timeVect(2) - timeVect(1))/2;
xVect = xVect(~idx);

idx = any(isinf(mean1) | isinf(mean2),2);
mean1 = mean1(~idx,:);
mean2 = mean2(~idx,:);

cDark = [0.2 0.1 0.6];
cLight = [0.8 0.9 0.8];



hPlot = plot(xVect,[nanmean(mean1);nanmean(mean2)]');
set(hPlot(1),'Color','black','LineWidth',1)
set(hPlot(2),'Color','black','LineWidth',1)
err1 = nanstd(mean1)./sqrt(sum(~isnan(mean1)));
err2 = nanstd(mean2)./sqrt(sum(~isnan(mean2)));
hPatch2 = patch([xVect fliplr(xVect)],...
    [nanmean(mean2)+err2 fliplr(nanmean(mean2)-err2)],...
    cDark);
hPatch1 = patch([xVect fliplr(xVect)],...
    [nanmean(mean1)+err1 fliplr(nanmean(mean1)-err1)],...
    cLight);



set([hPatch1 hPatch2],'FaceAlpha',0.5)
set(gca,'XLim',[0 xVect(end)])
xlabel('Inter turn interval (sec)')
ylabel('Switchiness')
legend([hPatch2,hPatch1],'Dark','Light')
legend('boxoff')

% clf
% plot((mean1./max(mean1,[],2))')

shg
% tVect = timeVect(1:(end-1));
% dat = nanmean(mean1);
%plot(tVect(1:120),fit_logistic(tVect(1:120),dat(1:120)))
%%
% MAPPING: Emax = b(1),  EC50 = b(2)
% hill_fit = @(b,x)  b(1).*x./(b(2)+x);
% b0 = [1, 8];
% B = lsqcurvefit(hill_fit, b0, tVect, nanmean(mean1));
% 
% fit1 = fit(timeVect(1:(end-1))',nanmean(mean1)','exp1')
% 
% [a,b] = glmfit(tVect',nanmean(mean1)','binomial','logit')
% logitFit = glmval(a,tVect,'logit');
% 
% plot(tVect,logitFit)
% 
% shg
