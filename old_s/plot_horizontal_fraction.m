function plot_horizontal_fraction(iviv_cells,L23h,L4h,L5h)
%% Plot horizontal fraction per layer
frac_h=L23h(iviv_cells,:);
frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position')
hold on;text(2.5,0.6,'L2/3','FontSize',12);
frac_h=L4h(iviv_cells,:);hold on;

frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
subplot(1,3,2);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');%ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');%ylabel('Fraction of total input');
box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position');
hold on;text(2.5,0.6,'L4','FontSize',12);

frac_h=L5h(iviv_cells,:);

frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
subplot(1,3,3);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');%ylabel('Fraction of total input');
box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');%ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position')
hold on;text(2.5,0.6,'L5','FontSize',12);
end