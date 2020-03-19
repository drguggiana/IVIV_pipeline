function plot_morphologies_iviv(str,ids,x,y)
fig7= figure;set(fig7, 'Name', 'Morphology');set(fig7, 'Position', [200, 0, 800, 800]);set(gcf,'color','w');   
 for i=1:length(ids)  
     subplot(x,y,i) 
     
    tmp=str{1,ids(i)}

       m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(tmp{1,2},[0 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(tmp{1,3},[0 0 1],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'

 set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
     tmp={};
 end

end