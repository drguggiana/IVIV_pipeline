function [stats_g] = display_inputs(frac_v,frac_h,frac_diffv,frac_diffh,groups)
set(0, 'DefaultAxesFontName', 'Arial')
%for all cells 
if isempty(groups)==1
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 600, 520]);set(gcf,'color','w');
%EX and IN Vertical
subplot(2,3,1);hold on;
for i=1:size(frac_v,1)
exp=plot(frac_v(i,1:16)',1:16,'-r');set(gca,'Ydir','reverse');ylabel('Vertical input');box off;
exp.Color(4) = 0.2;
exp=plot(frac_v(i,17:end)'*-1,1:16,'-b');set(gca,'Ydir','reverse');ylabel('Vertical map position');box off;
exp.Color(4) = 0.2;
end
hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
hold on;line([-1 1], [8 8],'Color','k','LineStyle','--');hold on;line([-1 1], [11 11],'Color','k','LineStyle','--');
hold on;mexp=errorbar(nanmean(frac_v(:,1:16)),1:16,nanstd(frac_v(:,1:16))/sqrt(size(frac_v,1)),'horizontal','k');set(gca,'Ydir','reverse');
hold on;mexp=errorbar(nanmean(frac_v(:,17:end))*-1,1:16,nanstd(frac_v(:,17:end))/sqrt(size(frac_v,1)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;xlim([-1 1]);xticks([-1:0.5:1]);hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});
hold on;text(-0.95,4.5,'L2/3');hold on;text(-0.95,7,'L4');hold on;text(-0.95,9.5,'L5');set(gca,'FontSize',10);
hold on;text(-0.5,0.2,'IN','Color','b');hold on;text(0.4,0.2,'EX','Color','r');
%EX and IN Vertical DIFFERENCE
subplot(2,3,2);hold on;
 for i=1:size(frac_diffv,1)
 exp=plot(frac_diffv(i,1:16)',1:16,'Color',[0.5 0.5 0.5]);set(gca,'Ydir','reverse');xlabel('Fraction of total input');box off;
 exp.Color(4) = 0.2;
 end
hold on;mexp=errorbar(nanmean(frac_diffv(:,1:16)),1:16,nanstd(frac_diffv(:,1:16))/sqrt(size(frac_diffv,1)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;xlim([-0.6 0.6]);
xticks([-0.6:0.3:0.6]);
hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});
hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
hold on;line([-1 1], [8 8],'Color','k','LineStyle','--');hold on;line([-1 1], [11 11],'Color','k','LineStyle','--');
set(gca,'FontSize',10);text(-0.3,0.2,'Difference','Color','k');

subplot(2,3,3);hold on;
mexp=errorbar(nanmean(frac_diffv(:,1:16))',1:16,nanstd(frac_diffv(:,1:16))/sqrt(size(frac_diffv,1)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;xlim([-0.2 0.2]);
xticks([-0.2:0.1:0.2]);hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});
hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
hold on;line([-1 1], [8 8],'Color','k','LineStyle','--');hold on;line([-1 1], [11 11],'Color','k','LineStyle','--');
set(gca,'FontSize',10);text(-0.1,0.2,'Difference','Color','k');


%EX and IN Horizontal
subplot(2,3,4);hold on;
for i=1:size(frac_v,1)
exp=plot(frac_h(i,1:16)','-r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
exp=plot(frac_h(i,17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-1 1]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');text(2.5,1,'Medial');text(10,1,'Lateral');set(gca,'FontSize',10);

%EX and IN Horizontal DIFFERENCE
subplot(2,3,5);hold on;
for i=1:size(frac_diffh,1)
exp=plot(frac_diffh(i,1:16),'Color',[0.5 0.5 0.5]);xlabel('Horizontal map position');box off;
exp.Color(4) = 0.2;
end
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'horizontal','k');
mexp.CapSize=3;
ylim([-0.4 0.4]);
yticks([-0.4:0.2:0.4]);
hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
stats_g=[];
hold on;text(5,0.42,'Difference','Color','k');

subplot(2,3,6);hold on;
;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'horizontal','k');
mexp.CapSize=3;
ylim([-0.1 0.1]);
yticks([-0.2:0.1:0.2]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
hold on;text(5,0.105,'Difference','Color','k');



%for grouping
else    
gr_nr=max(groups);
color_id={'r','g','b','k'}
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 400, 520]);set(gcf,'color','w');
%EX and IN Vertical
subplot(2,2,1);hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_v(find(groups==i),1:16))*-1,1:16,nanstd(frac_v(find(groups==i),1:16))/sqrt(size(frac_v(find(groups==i),1:16),1)),'horizontal',color_id{i});set(gca,'Ydir','reverse');
end
hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
hold on;line([-1 1], [9 9],'Color','k','LineStyle','--');xlim([-0.3 0.3]);
hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_v(find(groups==i),17:end)),1:16,nanstd(frac_v(find(groups==i),17:end))/sqrt(size(frac_v(find(groups==i),17:end),1)),'horizontal',color_id{i});set(gca,'Ydir','reverse');
end
%EX and IN Vertical difference
subplot(2,2,3);hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_diffv(find(groups==i),1:16))*-1,1:16,nanstd(frac_diffv(find(groups==i),1:16))/sqrt(size(frac_diffv(find(groups==i),1:16),1)),'horizontal',color_id{i});set(gca,'Ydir','reverse');
end
hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
hold on;line([-1 1], [9 9],'Color','k','LineStyle','--');xlim([-0.3 0.3]);

%EX and IN Horizontal
subplot(2,2,2);hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_h(find(groups==i),1:16))*-1,nanstd(frac_h(find(groups==i),1:16))/sqrt(size(frac_h(find(groups==i),1:16),1)),'horizontal',color_id{i});
end
hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_h(find(groups==i),17:end)),nanstd(frac_h(find(groups==i),17:end))/sqrt(size(frac_h(find(groups==i),17:end),1)),'horizontal',color_id{i});
end
mexp.CapSize=3;ylim([-0.3 0.3]);yticks([-0.3:0.1:0.3]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');

%EX and IN Horizontal difference
subplot(2,2,4);hold on;
for i=1:gr_nr
mexp=errorbar(nanmean(frac_diffh(find(groups==i),1:16))*-1,nanstd(frac_diffh(find(groups==i),1:16))/sqrt(size(frac_diffh(find(groups==i),1:16),1)),'horizontal',color_id{i});
hold on;
end
mexp.CapSize=3;ylim([-0.3 0.3]);yticks([-0.3:0.1:0.3]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');


%statistics
figure;
comb=[frac_v(:,3:10) groups'];
for i=1:size(comb,2)-1   
[p,tbl,stats] = kruskalwallis(comb(:,i),comb(:,end),'off');
multicom(:,:,i) = multcompare(stats,'CType','bonferroni');
signi_ev(i)=multicom(:,6,i)<0.05;
end

comb=[];
comb=[frac_v(:,19:26) groups'];
for i=1:size(comb,2)-1   
[p,tbl,stats] = kruskalwallis(comb(:,i),comb(:,end),'off');
multicom(:,:,i) = multcompare(stats,'CType','bonferroni');
signi_iv(i)=multicom(:,6,i)<0.05;
end

comb=[]
comb=[frac_h(:,5:12) groups'];
for i=1:size(comb,2)-1   
[p,tbl,stats] = kruskalwallis(comb(:,i),comb(:,end),'off');
multicom(:,:,i) = multcompare(stats,'CType','bonferroni');
signi_eh(i)=multicom(:,6,i)<0.05;
end

comb=[]
comb=[frac_h(:,21:28) groups'];
for i=1:size(comb,2)-1   
[p,tbl,stats] = kruskalwallis(comb(:,i),comb(:,end),'off');
multicom(:,:,i) = multcompare(stats,'CType','bonferroni');
signi_ih(i)=multicom(:,6,i)<0.05;
end

stats_g=[signi_ev;signi_iv;signi_eh;signi_ih]'
end