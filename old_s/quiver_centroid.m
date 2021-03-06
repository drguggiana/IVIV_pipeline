function quiver_centroid(ang1,ang2,ang3,idxtp,ex_map,in_map,t,j)
%t is the cell to plot as example
if j==1
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 1000, 200]);
%single example with vector
subplot(1,4,1);
exc_map=ex_map(:,:,t);
inplot_type = 2;sf=1;map_plot3(exc_map,'',inplot_type,fig1,sf,0,1);
hold on;plot(ang1(t,1),ang1(t,2),'k^','MarkerFaceColor','w','MarkerSize',5);
hold on;q=quiver(ang1(t,1),ang1(t,2),(ang1(t,3)-ang1(t,1)),(ang1(t,4)-ang1(t,2)),0);
q.Color='black';q.MaxHeadSize=0.9;q.LineWidth=1;
subplot(1,4,2);
inh_map=in_map(:,:,t);
inplot_type = 3;sf=1;map_plot3(inh_map,'',inplot_type,fig1,sf,0,1);
hold on;plot(ang2(t,1),ang2(t,2),'k^','MarkerFaceColor','w','MarkerSize',5);
hold on;q=quiver(ang2(t,1),ang2(t,2),(ang2(t,3)-ang2(t,1)),(ang2(t,4)-ang2(t,2)),0);
q.Color='black';q.MaxHeadSize=0.9;q.LineWidth=1;

subplot(1,4,3);
for i=1:length(idxtp)
%q=quiver(ang1(i,1),ang1(i,2),(ang1(i,3)-ang1(i,1)),(ang1(i,4)-ang1(i,2)),0);
q=quiver(8.5,ang1(i,2),(ang1(i,3)-ang1(i,1)),(ang1(i,4)-ang1(i,2)),0);
 q.Color='red';q.MaxHeadSize=0.5;q.LineWidth=1;
hold on;plot(8.5,ang1(i,2),'k^','MarkerFaceColor','w','MarkerSize',4);xlim([1 8]);ylim([1 8]);
%hold on;plot(c_cord(t,1),4,'k^');
end
for i=1:length(idxtp)
 hold on;plot(8.5,ang1(i,2),'^','MarkerFaceColor','w','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',4);
end
box off;set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([6 11]);ylim([1 8]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8.5 8.5],'Color','k','LineStyle','--');hold on;line([8.5 8.5], [1 16],'Color','k','LineStyle','--');axis off;

subplot(1,4,4);
for i=1:length(idxtp)
%q=quiver(ang2(i,1),ang2(i,2),(ang2(i,3)-ang2(i,1)),(ang2(i,4)-ang2(i,2)),0);
q=quiver(8.5,ang2(i,2),(ang2(i,3)-ang2(i,1)),(ang2(i,4)-ang2(i,2)),0);
 q.Color='blue';q.MaxHeadSize=0.5;q.LineWidth=1;
hold on;plot(8.5,ang2(i,2),'^','MarkerFaceColor','w','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',4);
%hold on;plot(c_cord(t,1),4,'k^');

end
for i=1:length(idxtp)
 hold on;plot(8.5,ang2(i,2),'^','MarkerFaceColor','w','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',4);
end
box off;set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([6 11]);ylim([1 8]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8.5 8.5],'Color','k','LineStyle','--');hold on;line([8.5 8.5], [1 16],'Color','k','LineStyle','--');axis off;

else
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]); 
 ran=[30 40 50 60 70 80 90 100 110 120 125 132];
 for i=1:12
     subplot(3,4,i)
hold on;line([8.5 8.5], [1 16],'Color','k','LineStyle','--');axis off;
 plot([8.5 ang1(ran(i),3)],[3.5 ang1(ran(i),4)],'r')
 hold on;plot(ang1(ran(i),3),ang1(ran(i),4),'ro');
 hold on;plot(8.5,3.5,'^','MarkerFaceColor','[0.5 0.5 0.5]','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',7);
 
 hold on;
 plot([8.5 ang2(ran(i),3)],[3.5 ang2(ran(i),4)],'b')
 hold on;plot(ang2(ran(i),3),ang2(ran(i),4),'bo');
  
 hold on;
 plot([8.5 ang3(ran(i),3)],[3.5 ang3(ran(i),4)],'r')
 hold on;scatter(ang3(ran(i),3),ang3(ran(i),4),'ro','filled');

box off;set(gca,'Ydir','reverse');

xlim([6 11]);ylim([1 8]);
hold on;
title(['#' num2str(ran(i))]);
 end
end
end