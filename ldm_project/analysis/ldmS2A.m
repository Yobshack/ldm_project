function ldmS2A(matFile_shi, matFile_trp)
%%
clf

matFile_shi = LDM_FigureS2AB_Shi;
matFile_trp = LDM_FigureS2AB_Trp;
hold on
set(gcf,'Renderer','painters')
marker_size = 5;
marker_alpha = 0.35;
pba = [12,1,1];
bg_color = [231,242,247]./255;

X_patch = (0:size(matFile_shi{:,1},1));
X_patch = interleave2(X_patch,X_patch);
ylimits = get(gca,'YLim');
Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
Y_patch = Y_patch(1:length(X_patch));



for ii = 1:length(genonames)
    % EXPERIMENTAL
    % SHIBIRE
    try
    shi_indx = contains(matFile_shi{4},matFile_shi{4}{ii});
    behave = matFile_shi{1}{shi_indx};
    idx = all(behave(:,1:2) > 0,2);
    activity = behave(idx,1:2);
    bias = behave(idx,3:4);
    labels = matFile_shi{4};
    
    jitter = rand(size(bias,1),1)*0.2;
    
    plt_x = (ii-1)+0.3 + jitter - mean(jitter);
    x_center = [ii-1+0.15,ii-1+0.45];
    color = [0 0 0];
    
    subplot(8,1,1)  
    pbaspect(pba) 
    hold on
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])

    
    scatter(plt_x,activity(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)  
    line(x_center,[median(activity(:,1)),median(activity(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
%     if ii == size(matFile_shi{:,1},1)
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
    
    subplot(8,1,2)
    hold on
    pbaspect(pba)
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
   

    
    scatter(plt_x,activity(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)
    line(x_center,[median(activity(:,2)),median(activity(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
    
    subplot(8,1,5)
    hold on
    pbaspect(pba)

    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,bias(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,1)),median(bias(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
%     
    subplot(8,1,6)
    hold on
    pbaspect(pba)

    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,bias(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,2)),median(bias(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
    
    % CONTROL
    behave = matFile_shi{2}{ii};
    idx = all(behave(:,1:2) > 0,2);
    activity = behave(idx,1:2);
    bias = behave(idx,3:4);
    
    jitter = rand(size(bias,1),1)*0.2;
    
    plt_x = (ii-1)-0.3 + 1  + jitter - mean(jitter);
    x_center = [ii-1-0.45+1,ii-0.15];
    color = [0.5 0.5 0.5];
    
    subplot(8,1,1)
    hold on
    pbaspect(pba)    

    scatter(plt_x,activity(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)  
    line(x_center,[median(activity(:,1)),median(activity(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
    subplot(8,1,2)
    hold on
    pbaspect(pba)
    scatter(plt_x,activity(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)
    line(x_center,[median(activity(:,2)),median(activity(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
    subplot(8,1,5)
    hold on
    pbaspect(pba)
    scatter(plt_x,bias(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,1)),median(bias(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
    subplot(8,1,6)
    hold on
    pbaspect(pba)

    scatter(plt_x,bias(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,2)),median(bias(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
    end
    
    % dTRPA1
    try
    trp_indx = contains(matFile_trp{4},matFile_shi{4}{ii});
    behave = matFile_trp{1}{trp_indx};
    idx = all(behave(:,1:2) > 0,2);
    activity = behave(idx,1:2);
    bias = behave(idx,3:4);
    labels = matFile_trp{4}{trp_indx};
    
    jitter = rand(size(bias,1),1)*0.2;
    
    plt_x = (ii-1)+0.3 + jitter - mean(jitter);
    x_center = [ii-1+0.15,ii-1+0.45];
    color = [0 0 0];
    
    subplot(8,1,3)  
    pbaspect(pba) 
    hold on
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,activity(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)  
    line(x_center,[median(activity(:,1)),median(activity(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
%     
    subplot(8,1,4)
    hold on
    pbaspect(pba)
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,activity(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)
    line(x_center,[median(activity(:,2)),median(activity(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
%     
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
    
    subplot(8,1,7)
    hold on
    pbaspect(pba)
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,bias(:,1),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,1)),median(bias(:,1))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',bg_color,'FaceAlpha',0.1,'EdgeColor','none')
%     end
    
    subplot(8,1,8)
    hold on
    pbaspect(pba)
    set(gca,'YMinorTick','off','TickLength',[0 0],'XTickLabelRotation',45,'XTick',[])
    
    scatter(plt_x,bias(:,2),'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',marker_alpha,'SizeData',marker_size)    
    line(x_center,[median(bias(:,2)),median(bias(:,2))],'LineWidth',3,...
        'Color',[1 0.5 0])
    
%     if ii == 52
%         ylimits = get(gca,'YLim');
%         Y_patch = repmat([ylimits([2,1]),ylimits([1,2])]',ceil(length(X_patch)/4),1);
%         Y_patch = Y_patch(1:length(X_patch));
%         area(X_patch',Y_patch','FaceColor',[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     end
    end
    
end


%print('FigS2.pdf','-dpdf','-bestfit')

dat = matFile_shi{4};
save('labelsS2AB.mat','dat');


%% Extracting stats for text

x = subplot(8,1,1);
y = subplot(8,1,3);
act_scores = [arrayfun(@(x) x.YData(1),x.Children(1:2:end)); arrayfun(@(x) x.YData(1),y.Children(1:2:end))];
avg_act = mean(act_scores);
min_act = min(act_scores);
max_act = max(act_scores);

x = subplot(8,1,2);
y = subplot(8,1,4);
act_scores = [arrayfun(@(x) x.YData(1),x.Children(1:2:end)); arrayfun(@(x) x.YData(1),y.Children(1:2:end))];
avg_act = mean(act_scores);
min_act = min(act_scores);
max_act = max(act_scores);

x = subplot(8,1,5);
y = subplot(8,1,7);
tb_scores = [arrayfun(@(x) nanmean(x.YData),x.Children(2:2:end));...
    arrayfun(@(x) nanmean(x.YData),y.Children(2:2:end))];
avg_sd = mean(tb_scores);
min_sd = min(tb_scores);
max_sd = max(tb_scores);

x = subplot(8,1,6);
y = subplot(8,1,8);
act_scores = [arrayfun(@(x) x.YData(1),x.Children(1:2:end)); arrayfun(@(x) x.YData(1),y.Children(1:2:end))];
avg_act = mean(act_scores);
min_act = min(act_scores);
max_act = max(act_scores);


%%