function [statsout]=extent_diff(diffL23fr,diffL4fr,diffL5fr,span);
%histograms for difference in L23, L4, L5 Vertical and Horizontal (span),
%outputs ttests for comparing across and within layers with

ex_spanhL23=span(:,1);
ex_spanhL4=span(:,2);
ex_spanhL5=span(:,3);
in_spanhL23=span(:,4);
in_spanhL4=span(:,5);
in_spanhL5=span(:,6);

xf=69;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 400]);
subplot(3,2,1);
;hold on;
h4 = histogram(diffL23fr,8);h4.BinWidth = 0.07;h4.EdgeColor = 'k';h4.FaceColor = 'm';hold on;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylim([0 50]);yticks([0:25:50])
title('Vertical');ylabel('Cell counts');hold on;set(gca,'FontSize',10);
xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);
text(-0.55,50,'IN','Color','b');text(0.55,50,'EX','Color','r');text(-0.55,40,'L2/3','Color','k');
hold on;plot(nanmean(diffL23fr),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');
subplot(3,2,3);
h4 = histogram(diffL4fr,8);h4.BinWidth =  0.07;h4.EdgeColor = 'k';h4.FaceColor = 'g'
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylabel('Cell counts');hold on;set(gca,'FontSize',10);box off;
text(-0.55,40,'L4');
hold on;plot(nanmean(diffL4fr),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);
subplot(3,2,5);
h4 = histogram(nonzeros(diffL5fr),8);h4.BinWidth =  0.07;h4.FaceColor = [0.5 0.5 0.5];ylabel('Cell counts');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');
xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);xlabel('\Delta EX-IN vertical fraction');hold on;set(gca,'FontSize',10);
text(-0.55,40,'L5','Color','k');
hold on;plot(nanmean(diffL5fr),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
%Horizontal span
subplot(3,2,2);
h4 = histogram((ex_spanhL23-in_spanhL23)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'm';
xlim([-10*xf 10*xf]);hold on; title('Horizontal');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,50,'IN','Color','b');text(7*xf,50,'EX','Color','r');
text(-7*xf,40,'L2/3','Color','k');xlim([-600 600]);xticks([-600:300:600]);
hold on;plot(nanmean((ex_spanhL23-in_spanhL23)*xf),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,4);
hold on;h4 = histogram((ex_spanhL4-in_spanhL4)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'g';
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,40,'L4','Color','k');xlim([-600 600]);xticks([-600:300:600]);
hold on;plot(nanmean((ex_spanhL4-in_spanhL4)*xf),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,6);
hold on;h4 = histogram((ex_spanhL5-in_spanhL5)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];;
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
xlabel('\Delta EX-IN horizontal span (µm)');hold on;set(gca,'FontSize',10);
text(-7*xf,40,'L5','Color','k');xlim([-600 600]);xticks([-600:300:600]);
hold on;plot(nanmean((ex_spanhL5-in_spanhL5)*xf),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])

%Statistics
temp=[diffL23fr diffL4fr diffL5fr ex_spanhL23-in_spanhL23 ex_spanhL4-in_spanhL4 ex_spanhL5-in_spanhL5];
for i=1:6
[r(i) k(i)]=ttest(temp(:,i));
end

[p1,tbl,stats1] = anova1(temp(:,1:3),[],'off');
multicom = multcompare(stats1,'CType','bonferroni','Display','off');
[p1,tbl,stats2] = anova1(temp(:,4:6),[],'off');
multicom2 = multcompare(stats2,'CType','bonferroni','Display','off');

statsout.ttests=k;
statsout.vert=multicom(:,6)';
statsout.span=multicom2(:,6)';
end