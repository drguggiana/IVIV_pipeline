function shift_ori_oe(fit_ex1,fit_ex1oe,fit_ex2,fit_ex2oe,ind_tr1,ind_tr1oe,ind_tr2,ind_tr2oe,ori1,ori2,tw1,tw2)
ori=90-[1:1:180];
pori=[0:45:315];
peak=max(fit_ex1);
find(fit_ex1==peak);
peak2=max(fit_ex2);
find(fit_ex2==peak2);

format short g;

if isempty(ind_tr1)==0
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 200, 700, 250]);
subplot(1,3,1)
%errorbar(nanmean(ind_tr1,2)/max(nanmean(ind_tr1,2)),(nanstd(ind_tr1,[],2)/sqrt(4))/max(nanmean(ind_tr1,2)))
errorbar(pori,nanmean(ind_tr1,2)/max(nanmean(ind_tr1,2)),nanstd(ind_tr1/max(nanmean(ind_tr1,2)),[],2)/sqrt(4),'-o','Color','b');
hold on;
errorbar(pori,nanmean(ind_tr1oe,2)/max(nanmean(ind_tr1,2)),nanstd(ind_tr1oe/max(nanmean(ind_tr1,2)),[],2)/sqrt(4),'-o','Color','r');
% plot(fit_ex1/max(fit_ex1),'k');
hold on;
text(140,1.4,['OSI= ' num2str(round(ori1,2))],'Color','b');
hold on;box off
ylabel('Normalized response');xlabel('Orientation (deg)');
yticks([0:0.5:1.5]);
ylim([0 1.5]);
xlim([0 315]);xticks([0:45:315]);xtickangle(45);
set(gca,'FontSize',10);

subplot(1,3,2)
%plot(fit_ex2/max(fit_ex2),'b');box off;
errorbar(pori,nanmean(ind_tr2,2)/max(nanmean(ind_tr2,2)),nanstd(ind_tr2/max(nanmean(ind_tr2,2)),[],2)/sqrt(4),'-o','Color','r');
hold on; text(140,1.4,['OSI= ' num2str(round(ori2,2))],'Color','r');
hold on;
errorbar(pori,nanmean(ind_tr2oe,2)/max(nanmean(ind_tr2,2)),nanstd(ind_tr2oe/max(nanmean(ind_tr2,2)),[],2)/sqrt(4),'-o','Color','b');
ylabel('Normalized response');xlabel('Orientation (deg)');
ylim([0 1.5]);yticks([0:0.5:1.5]);
xlim([0 315]);xticks([0:45:315]);xtickangle(45);
set(gca,'FontSize',10);
box off;

subplot(1,3,3)
plot(circshift(ori,find(fit_ex1==peak)-90),fit_ex1/max(fit_ex1),'k');
hold on;
plot(circshift(ori,find(fit_ex2==peak2)-90),fit_ex2/max(fit_ex2),'b');
hold on;box off;
xlabel('Distance from pref. orientation (deg)');
text(50,0.3,['TW= ' num2str(round(tw1,2))]);
hold on;box off
text(50,0.2,['TW= ' num2str(round(tw2,2))],'Color','b');
yticks([0:0.2:1]);
set(gca,'FontSize',10);

else
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 200, 250]);   
plot(circshift(ori,find(fit_ex1==peak)-90),fit_ex1/max(fit_ex1),'k');
hold on;
plot(circshift(ori,find(fit_ex2==peak2)-90),fit_ex2/max(fit_ex2),'b');
hold on;box off;
xlabel('Distance from preferred orientation (deg)');
text(50,0.25,['TW= ' num2str(tw1)]);
hold on;box off
text(50,0.2,['TW= ' num2str(tw2)],'Color','b');
yticks([0:0.2:1]);
set(gca,'FontSize',12);
end
end