%% Kyobi Skutt-Kakaria
% 02.04.2018 - Created
% 02.04.2018 - updated

% this function generates kdes for low to high temp shifts for each of the 5 primary metric
% differences for the paper

function out = generateShiftPlot(array,data1,error1,behave)



YNames = {'Activity','Turn Bias', 'Mutual Information','Wall Distance','Inter-turn Interval Variability'};
XNames = {'turns per minute','proportion right Turns', 'normalized probability right follows right',...
    'wall proximity','inter-turn interval variability'};

for ii = behave
    %subplot(2,5,ii)
    hold on
    if size(data1,2) == 20
    dataHigh = [data1(:,ii+10), data1(:,ii+15)];
    dataLow = [data1(:,ii), data1(:,ii+5)];
    else
    dataHigh = [data1(:,1:2)];
    dataLow = [data1(:,3:4)];
    end
    idx = any(data1(:,[1 6]) < 0.03,2);
    
    dataHigh(idx,:) = nan;
    dataLow(idx,:) = nan;
    
    f = @(x) (x - nanmean(x,1))./nanstd(x,[],1);
    dataHigh = f(dataHigh);
    dataLow = f(dataLow);
    
    idxISO = strcmp(array,'ISO');
    clear lightISO lightSHI
    if sum(idxISO) > 0
        dataISO = [dataHigh(idxISO,:), dataLow(idxISO,:)];
        dataISO = dataISO(all(~isnan(dataISO),2),:);
        
        [~,temp] = min(abs(dataISO(:,1) - dataISO(:,3:4)),[],2);
        lightISO(:,1) = temp == 1;
        
        [~,temp] =  min(abs(dataISO(:,2) - dataISO(:,3:4)),[],2);
        lightISO(:,2) = temp == 1;
        
        lightISO = double(lightISO);
        
        lightISO(lightISO == 0) = -1;
        
        proportion(1,:) = nanmean(lightISO,1);
        
        for jj = 1:1000
            idx = randi(size(dataISO,1),size(dataISO,1),1);
%             [a,b] = sort(dataISO(idx,:));
%             [b1, i1] = sort(b);
%             [~,temp] = min(abs(i1(:,1) - i1(:,3:4)),[],2);
%             lightISO(:,1) = temp == 1;
%             
%             [~,temp] = min(abs(i1(:,2) - i1(:,3:4)),[],2);
%             lightISO(:,2) = temp == 1;
% % 
            [~,temp] = min(abs(dataISO(idx,1) - dataISO(idx,3:4)),[],2);
            lightISO(:,1) = temp == 1;
            
            [~,temp] =  min(abs(dataISO(idx,2) - dataISO(idx,3:4)),[],2);
            lightISO(:,2) = temp == 1;
            
            lightISO = double(lightISO);
            
            lightISO(lightISO == 0) = -1;
            
            reSamp1(jj,:) = nanmean(lightISO,1);
        end
        proportion(1,:) = nanmean(reSamp1,1);
        proportionE(1,:) = nanstd(reSamp1,[],1);
    end
    
    idxSHI = strcmp(array,'SHI');
        if sum(idxSHI) > 0


        dataSHI = [dataHigh(idxSHI,:), dataLow(idxSHI,:)];
        dataSHI = dataSHI(all(~isnan(dataSHI),2),:);
        
        [~,temp] = min(abs(dataSHI(:,1) - dataSHI(:,3:4)),[],2);
        lightSHI(:,1) = temp == 1;
        [~,temp] =  min(abs(dataSHI(:,2) - dataSHI(:,3:4)),[],2);
        lightSHI(:,2) = temp == 1;
        
        lightSHI = double(lightSHI);
        lightSHI(lightSHI == 0) = -1;
        proportion(2,:) = nanmean(lightSHI,1);
        

        
        for jj = 1:1000
            idx = randi(size(dataSHI,1),size(dataSHI,1),1);
%             [a,b] = sort(dataSHI(idx,:));
%             [b1, i1] = sort(b);
%             [~,temp] = min(abs(i1(:,1) - i1(:,3:4)),[],2);
%             lightSHI(:,1) = temp == 1;
%             
%             [~,temp] = min(abs(i1(:,2) - i1(:,3:4)),[],2);
%             lightSHI(:,2) = temp == 1;
% %             
            [~,temp] = min(abs(dataSHI(idx,1) - dataSHI(idx,3:4)),[],2);
            lightSHI(:,1) = temp == 1;
            
            [~,temp] =  min(abs(dataSHI(idx,2) - dataSHI(idx,3:4)),[],2);
            lightSHI(:,2) = temp == 1;
% %             
            lightSHI = double(lightSHI);
            
            lightSHI(lightSHI == 0) = -1;
            
            reSamp2(jj,:) = nanmean(lightSHI,1);
        end
        proportionE(2,:) = nanstd(reSamp2,[],1);
    end
    
    %
    %     idxTRP = strcmp(array,'TRP');
    %     if sum(idxTRP) > 0
    %     dataTRP = [dataHigh(idxTRP,:), dataLow(idxTRP,:)];
    %
    %     diffs = abs(dataTRP(:,1) - dataTRP(:,2));
    %     absDev = abs(diffs - nanmean(diffs)) > nanstd(diffs);
    %     dataTRP = dataTRP(all(~isnan(dataTRP),2),:);
    %
    %     for jj = 1:1000
    %         idx = randi(size(dataTRP,1),size(dataTRP,1),1);
    %         [a,b] = sort(dataTRP(idx,:));
    %         [b1, i1] = sort(b);
    %         [~,temp] = min(abs(i1(:,1) - i1(:,3:4)),[],2);
    %         lightTRP(:,1) = temp == 1;
    %
    %         [~,temp] = min(abs(i1(:,2) - i1(:,3:4)),[],2);
    %         lightTRP(:,2) = temp == 1;
    %
    %         reSamp3(jj,:) = sum(lightTRP)./size(lightTRP,1);
    %     end
    %     proportion(3,:) = nanmean(reSamp3);
    %     proportionE(3,:) = nanstd(reSamp3);
    %     end

    
    %%
    
    %     c1 = [0.4 0.3 0.7];
    %     c2 = [0.8 0.9 0.8];
    %
    %     clf
    %
    %     hBar = bar(flipud([dat 1-dat]),'BarLayout','stacked','FaceAlpha',1,...
    %         'BarWidth',1);
    %     hBar(1).FaceColor = c2;
    %     hBar(2).FaceColor = c1;
    %     set(gca,'YTickLabel',{'TRP : Dark', 'TRP : Light','SHI : Dark',...
    %         'SHI : Light','ISO : Dark','ISO : Light'})
    %     line([0.5 0.5], [0 7],'LineWidth',2,'Color','red')
    %     hold on
    %     errorbar(flipud(dat),1:length(dat),flipud(err),...
    %         'horizontal','LineStyle','none','CapSize',0,'Color','black','LineWidth',5)
    %     shg
    %%
    
    
    
    y = proportion';
    e = proportionE';
    n = {num2str(size(dataISO,1)),num2str(size(dataSHI,1))};
    hold on

    hBar = bar(y);
    hBar(1).FaceColor = [0/255 114/255 189/255];
    hBar(2).FaceColor = [217/255 83/255 25/255];
    drawnow
    for ii = 1:size(y,1)
        hError = errorbar(hBar(ii).XData+hBar(ii).XOffset,y(:,ii),e(:,ii),...
            'LineStyle','none','Color','black','CapSize',0,'LineWidth',3);
    end
    
    shg
    
   
    if size(y,2) ==2
        for jj = 1:size(y,1)
            pvals{jj} = num2str(flyVacPersMetricPValue([y(jj,1),e(jj,1),y(jj,2),e(jj,2)]));
        end
        text(2,0.6,{pvals{:},n{:}})
    end
    
    %     x = normrnd(dat(1),err(1),size(dataISO,1),1);
    %     y = normrnd(dat(3),err(3),size(dataSHI,1),1);
    %     [h,p1] = ttest2(x,y,'Vartype','unequal');
    %     x = normrnd(dat(2),err(2),size(dataISO,1),1);
    %     y = normrnd(dat(4),err(4),size(dataSHI,1),1);
    %     [h,p1] = ttest2(x,y,'Vartype','unequal');
    %     x = normrnd(dat(1),err(1),size(dataISO,1),1);
    %     y = normrnd(dat(5),err(5),size(dataTRP,1),1);
    %     [h,p1] = ttest2(x,y,'Vartype','unequal');
    %         x = normrnd(dat(2),err(2),size(dataISO,1),1);
    %     y = normrnd(dat(6),err(6),size(dataTRP,1),1);
    %     [h,p1] = ttest2(x,y,'Vartype','unequal');
    %
    
    set(gca,'YLim',[-0.6 0.6],'XTick',[1:2],'XTickLabel',{'Light','Dark'},'FontSize',24,...
        'FontWeight','bold','FontName','arial','YTick',-1:0.2:1,...
        'Color','none','XLim',[0.5 2.5],'Box','on')
    set(gcf,'Renderer','painters')
    %colormap(brewermap(12,'paired'))
    pbaspect([1 2 1])
    
    line([0 3],[0 0],'Color',[0.7 0 0],'LineWidth',2)
    legend({'Gal4/+','Gal4/Shits'})
    legend('boxoff')
    % %%
    % clf
    % hold on
    %
    % data = [nanmean(resampISO) nanmean(resampSHI) nanmean(resampTRP)];
    % error = [nanstd(resampISO) nanstd(resampSHI) nanstd(resampTRP)];
    %
    % hBar1 = bar(data(1:2:6),'BarLayout','grouped');
    % hBar2 = bar(data(2:2:6),'BarLayout','grouped');
    %
    % set([hBar1 hBar2],'BarWidth',0.2,'LineWidth',1)
    % set(hBar1,'FaceColor',[0.2 0.1 0.6],'XData',[0.4 1.1 1.8],'FaceAlpha',0.8)
    % set(hBar2,'FaceColor',[0.8 0.9 0.8],'XData',[0.55 1.25 1.95],'FaceAlpha',0.8)
    %
    % hErr = errorbar(interleave2(hBar1.XData,hBar2.XData)',data,error,'LineStyle','none');
    % set(hErr,'CapSize',0,'Color','black','LineWidth',6)
    %
    % set(gca,'FontSize',20,'FontWeight','bold','XLim',[0 1.5],'YLim',[-0.5 0.5],'XTick',[],...
    %     'XLim',[0 2.5])
    % shg
    %
    %
    % %[~,p1] = ttest(shiftL,dataNull);
    % %[~,p2] = ztest(nanmean(resamp(:,1)),0,nanstd(resamp(:,1)));
    %
    % % vbes = [nanmean(IBE)];
    % % vbee = [nanmean(IBE)];
    %
    % % tStr = strcat({'med_{delta} = '}, num2str(round(nanmedian(data),2)),...
    % %     {'\newlineIBE_{delta} = '}, num2str(round(vbes(1),3)));
    % % hText = text(xMax*0.3, yMax*0.6,tStr,'Interpreter','tex','FontWeight','bold','FontSize',12);
    % %
    % % if p1 < 0.001
    % %     sym1 = '***';
    % % elseif p1 < 0.01
    % %     sym1 = '**';
    % % elseif p1 < 0.05
    % %     sym1 = '*';
    % % else
    % %     sym1 = '';
    % % end
    % %
    % % if p2 < 0.001
    % %     sym2 = '***';
    % % elseif p2 < 0.01
    % %     sym2 = '**';
    % % elseif p2 < 0.05
    % %     sym2 = '*';
    % % else
    % %     sym2 = '';
    % % end
    % % %
    % % % p1 = hText.Position(1)*2.2;
    % %
    % % text(0.25,0.3,sym2,'FontSize',14,...
    % %     'HorizontalAlignment','center','FontWeight','bold')
    % %
    % % text(p1,p2*0.94,sym2,'FontSize',14,...
    % %     'HorizontalAlignment','center','FontWeight','bold')
    %
    % pbaspect([1 1 1])clf
    %
    % hLegend = legend([hBar1,hBar2],{'Dark','Light'},'Location','northeast');
    % legend('boxoff')
    %
    % hTitle = title(YNames(ii));
    %
    % ylabel('\rho_{light} - \rho_{dark}','Rotation',90,'VerticalAlignment','cap')
    %
    % drawnow
    % shg
end

out = [y,e]

