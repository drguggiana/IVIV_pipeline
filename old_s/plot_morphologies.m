function plot_morphologies(str,ids,co_a,x,y)

 fig7= figure;set(fig7, 'Name', 'Morphology');set(fig7, 'Position', [200, 0, 400, 400]);set(gcf,'color','w');   
 for i=1:length(ids)  
     subplot(x,y,i)
     
     if ~isempty(str(ids(i)).morphtraces)==1
         if co_a==1
       m=plot_tree(str(ids(i)).morphtraces{1,1},[1 0 0],[0 str(ids(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'k';
         else
        m=plot_tree(str(ids(i)).morphtraces{1,1},[1 0 0],[0 str(ids(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'b';
         end
       m1=plot_tree(str(ids(i)).morphtraces{1,2},[0 0 0],[0 str(ids(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = [0.5 0.5 0.5];
       m2=plot_tree(str(ids(i)).morphtraces{1,3},[0 0 1],[0 str(ids(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = [0.5 0.5 0.5];
else
    m=plot(str(ids(i)).pialD,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
end
 set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
     
 end


end