function plot_morphologies_individual(str,ids,only_apical,co_a,x,y)
%figure;
 fig7= figure;set(fig7, 'Name', 'Morphology');set(fig7, 'Position', [200, 0, 600, 800]);set(gcf,'color','w');   
 for i=1:length(ids)  
     subplot(x,y,i)
     
     if ~isempty(str(ids(i)).tr_apical)==1
         if only_apical==1
         if co_a==1
       m=plot_tree(str(ids(i)).tr_apical,[1 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
        m.EdgeColor = 'k';
         else
        m=plot_tree(str(ids(i)).tr_apical,[1 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
        m.EdgeColor = 'b';
         end
         elseif only_apical==0
       m1=plot_tree(str(ids(i)).tr_basal,[0 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
       m1.EdgeColor = [0.5 0.5 0.5];
       m1.EdgeColor = 'm';
       m2=plot_tree(str(ids(i)).tr_soma,[0 0 1],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
       m2.EdgeColor = [0.5 0.5 0.5];
         else only_apical==2
              if co_a==1
       m=plot_tree(str(ids(i)).tr_apical,[1 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
        m.EdgeColor = 'k';
               else
        m=plot_tree(str(ids(i)).tr_apical,[1 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
        m.EdgeColor = 'b';
              end
        m1=plot_tree(str(ids(i)).tr_basal,[0 0 0],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
       m1.EdgeColor = [0.5 0.5 0.5];
       m1.EdgeColor = 'm';
       m2=plot_tree(str(ids(i)).tr_soma,[0 0 1],[0 str(ids(i)).pia 0],[],1,'-b');hold on;
       m2.EdgeColor = [0.5 0.5 0.5];
         
         end
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