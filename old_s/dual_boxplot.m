function [statsout]=dual_boxplot(par,g1,g2,t)
gr_nr=2;
%color_id={[0.7 0 0.4],[0 0.5 0.5]};
%color_id={[0.5 0.5 0.5],[0.7 0 0.4]};
%color_id={[0.5 0.5 0.5],[0.8500, 0.3250, 0.0980]};
%color_id={[0 0 1],[1 0 0]};
%color_id={[1 0 0],[0.5 0.5 0.5]};
color_id={[0 0 0],[0.8500 0.3250 0.0980]};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 200]);set(gcf,'color','w');
for i=1:gr_nr
gr_m=[nanmean(par(g1)) nanmean(par(g2))];  

gr_sem=[nanstd(par(g1))/sqrt(sum(~isnan(par(g1)))) nanstd(par(g2))/sqrt(sum(~isnan(par(g2))))];
hold on;
plot(ones(1,length(par(g1))),par(g1),'o','MarkerEdgeColor',color_id{1},'MarkerSize',5);
hold on;
plot(ones(1,length(par(g2)))*2,par(g2),'o','MarkerEdgeColor',color_id{2},'MarkerSize',5);


if t==0

[p k]=ranksum(par(g1),par(g2));
    statsout=p;
elseif t==1
    [p k]=ranksum(par(g1),par(g2),'tail','right');
    statsout=p;
elseif t==2
    [p k]=ranksum(par(g1),par(g2),'tail','left');
    statsout=p;
else    [k p]=ttest2(par(g1),par(g2));
    statsout=p;
end

hold on;
title(['p=' num2str(p)],'FontWeight','normal');
end
for i=1:gr_nr
hold on;
er=errorbar(i,gr_m(i),gr_sem(i));er.Color = color_id{i};er.LineWidth=1.5;er.LineStyle = 'none'; 
er.CapSize = 10;hold on;
hold on;plot(i,gr_m(i),'o','MarkerFaceColor',color_id{i},'MarkerEdgeColor','k','MarkerSize',7);
end
end