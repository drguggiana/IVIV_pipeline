function [statsout] = barplot_sw(data,groups_idx,labelxy)
gr_nr=max(groups_idx);

%color_id={[0.5 0.5 0.5],[0.7 0 0.4],[0 0.5 0.5]};
color_id={'w','w','w','w','w','w'};
%clus_id={'C1','C2','C3','C4','C5'};
%clus_id={'ipsi','contra','bino'};
%clus_id={'non resp.','10-60°','100-150°'};
clus_id={'250','10','1','250','10','1'};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [200, 0, 300, 300]);set(gcf,'color','w');


for i=1:gr_nr
gr_m(i)=nanmean(data(find(groups_idx==i)));  
gr_sem(i)=nanstd(data(find(groups_idx==i)))/sqrt(length(find(groups_idx==i)));
%plotting
hold on;
%bar plot
b=bar(i,gr_m(i));b.FaceColor=color_id{i};
%individual data points
hold on;
plot(ones(1,length(~isnan(data(find(groups_idx==i)))))*i,data(find(groups_idx==i)),'ko','MarkerEdgeColor',[0.7,0.7,0.7])
%errorbars
hold on;
er=errorbar(i,gr_m(i),gr_sem(i));er.Color = [0 0 0];er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
hold on;xticks([1:1:gr_nr]);hold on;set(gca,'XTickLabel',{clus_id{1:gr_nr}});hold on;
xlabel(labelxy{1});ylabel(labelxy{2});
%text(1:length(c_size),c_size,num2str(c_size'),'vert','bottom','horiz','center');
box off;
end
set(gca,'FontSize',10);
%-----------------------------------------------% statistics---------------
comb=[data groups_idx];
figure();
[p,tbl,stats] = kruskalwallis(comb(:,1),comb(:,end),'off');
multicom = multcompare(stats,'CType','bonferroni');
statsout=multicom(:,6);
end