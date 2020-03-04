function quiver_centroid(ang1,ang2,idxtp,ex_map,in_map,t)
%t is the cell to plot as example

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 200]);
%single example with vector
subplot(1,3,1);
exc_map=ex_map(:,:,t);
inplot_type = 2;sf=1;map_plot3(exc_map,'',inplot_type,fig1,sf,0,1);
hold on;plot(ang1(t,1),ang1(t,2),'k^','MarkerFaceColor','w','MarkerSize',3);
hold on;q=quiver(ang1(t,1),ang1(t,2),(ang1(t,3)-ang1(t,1)),(ang1(t,4)-ang1(t,2)),0);
q.Color='black';q.MaxHeadSize=0.8;q.LineWidth=1;
subplot(1,3,2);
inh_map=in_map(:,:,t);
inplot_type = 3;sf=1;map_plot3(inh_map,'',inplot_type,fig1,sf,0,1);
hold on;plot(ang2(t,1),ang2(t,2),'k^','MarkerFaceColor','w','MarkerSize',3);
hold on;q=quiver(ang2(t,1),ang2(t,2),(ang2(t,3)-ang2(t,1)),(ang2(t,4)-ang2(t,2)),0);
q.Color='black';q.MaxHeadSize=0.8;q.LineWidth=1;
subplot(1,3,3);
for i=1:length(idxtp)
q=quiver(ang1(i,1),ang1(i,2),(ang1(i,3)-ang1(i,1)),(ang1(i,4)-ang1(i,2)),0);
 q.Color='red';q.MaxHeadSize=0.1;q.LineWidth=1;
hold on;plot(ang1(i,1),ang1(i,2),'k^','MarkerFaceColor','w','MarkerSize',3);
%hold on;plot(c_cord(t,1),4,'k^');
end
hold on;
for i=1:length(idxtp)
q=quiver(ang2(i,1),ang2(i,2),(ang2(i,3)-ang2(i,1)),(ang2(i,4)-ang2(i,2)),0);
 q.Color='blue';q.MaxHeadSize=0.1;q.LineWidth=1;
hold on;plot(ang2(i,1),ang2(i,2),'k^','MarkerFaceColor','w','MarkerSize',3);
%hold on;plot(c_cord(t,1),4,'k^');
end
box off;set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([1 16]);ylim([1 16]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8 8],'Color','k','LineStyle','--');hold on;line([8 8], [1 16],'Color','k','LineStyle','--');axis off;

end