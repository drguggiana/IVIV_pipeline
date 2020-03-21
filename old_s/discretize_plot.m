function discretize_plot(dis_val,binnr,dvalue,avgt)

dis=discretize(horzcat(dis_val),binnr)
if avgt==1
for i=1:binnr
    idxgr=find(dis==i);
    g_avg(i)=nanmean(dvalue(idxgr));
    g_sem(i)=nanstd(dvalue(idxgr))/sqrt(length(idxgr));
    idxgr=[];
end
figure;set(gcf,'color','w');
errorbar(g_avg,g_sem,'-o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');box off;
else avgt==0

end
end