function plot_morphologies_iviv(str,ids,x,y,m_flip_a)
fig7= figure;set(fig7, 'Name', 'Morphology');set(fig7, 'Position', [200, 0, 800, 800]);set(gcf,'color','w');   
 for i=1:length(ids)  
     subplot(x,y,i) 
     
    tmp=str{1,ids(i)}
   
if iscell(tmp)==1
       m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'k'
       m1=plot_tree(tmp{1,2},[0 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'm'
       m2=plot_tree(tmp{1,3},[0 0 1],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'g'

 set(gca,'Ydir','reverse');
  xlim([-300 300]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
     tmp={};
     if m_flip_a(ids(i))==1
         set(gca,'Xdir','reverse');
     else
     end
else
     plot(1,1,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
      tmp={};axis off
end
 end

 figure;set(gcf,'color','w');
 hold on
 for i=1:length(ids)  
     tmp=str{1,ids(i)}
  if iscell(tmp)==1
       m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'k'
       m1=plot_tree(tmp{1,2},[0 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'm'
       m2=plot_tree(tmp{1,3},[0 0 1],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'g'

 set(gca,'Ydir','reverse');
  xlim([-300 300]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
     tmp={};
     hold on;
      if m_flip_a(ids(i))==1
         set(gca,'Xdir','reverse');
     else
     end
else
    
      tmp={};axis off
end   
 end
 
end