function [statsout]=dual_barplot(par,g1,g2,t)
gr_nr=2;
color_id={[0.7 0 0.4],[0 0.5 0.5]};
%color_id={[0.5 0.5 0.5],[0.7 0 0.4]};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 100, 200]);set(gcf,'color','w');
for i=1:gr_nr
gr_m=[nanmean(par(g1)) nanmean(par(g2))];  
gr_sem=[nanstd(par(g1))/sqrt(length(par(g1))) nanstd(par(g2))/sqrt(length(par(g2)))];
hold on;
b=bar(i,gr_m(i));b.FaceColor=color_id{i};
hold on;
plot(ones(1,length(par(g1))),par(g1),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
hold on;
plot(ones(1,length(par(g2)))*2,par(g2),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);


if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par(g1),par(g2));
    statsout=p;
elseif t==1
    [p k]=ranksum(par(g1),par(g2),'tail','right');
    statsout=p;
else t==2
    [p k]=ranksum(par(g1),par(g2),'tail','left');
    statsout=p;
end

hold on;
title(['p=' num2str(p)],'FontWeight','normal');
end
for i=1:gr_nr
hold on;
er=errorbar(i,gr_m(i),gr_sem(i));er.Color = [0 0 0];er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
end
end