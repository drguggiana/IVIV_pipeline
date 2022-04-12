function paired_plot_box(data,cl)



fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 150, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
%hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','r','Markersize',7);
hold on;er1=errorbar([0.75],nanmean(data(:,1)),nanstd(data(:,1),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor',cl{1},'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(data(:,2)),nanstd(data(:,2),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor',cl{2},'Markersize',8);
 [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');

end