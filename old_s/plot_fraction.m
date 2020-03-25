function plot_fraction(L23fr,L4fr,L5fr,pia_input)

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 450, 500]);
subplot(3,2,1);
par1=L23fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;%xlabel('Input fraction');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])  
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['L23: r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'k-');set(gca,'box','off');set(gcf,'color','w');;   set(gca,'Ydir','reverse');
    ylabel('Pial depth (µm)'); set(gca,'FontSize',10);%text(0.05,120,'L2/3');
    
subplot(3,2,3);
par1=L4fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['L4: r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'k-');set(gca,'box','off');set(gcf,'color','w');   set(gca,'Ydir','reverse');
     ylabel('Pial depth (µm)');%text(0.05,120,'L4');set(gca,'FontSize',10)
 set(gca,'FontSize',10)
 
 subplot(3,2,5);   
 par1=L5fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['L5: r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
 
    hold on;
    %plot(par1,yfit,'k-');set(gca,'box','off');set(gcf,'color','w');;  
    set(gca,'Ydir','reverse');
     ylabel('Pial depth (µm)');%text(0.05,110,'L5');set(gca,'FontSize',10);
     xlim([0. 0.4]);xlabel('Fraction total input');set(gca,'FontSize',10);

 subplot(3,2,2);
par1=L23fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;%xlabel('Input fraction');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])  
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['L23: r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end  
    set(gca,'Ydir','reverse');
  set(gca,'FontSize',10);%text(0.05,120,'L2/3');
 
 subplot(3,2,4);
par1=L4fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;%xlabel('Input fraction');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])  
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['L4: r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end  
 P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'k-');set(gca,'box','off');set(gcf,'color','w');  set(gca,'Ydir','reverse');
     set(gca,'Ydir','reverse');
set(gca,'FontSize',10);  xlim([0. 0.4]);%text(0.05,120,'L4'); 


 subplot(3,2,6);
par1=L5fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;%xlabel('Input fraction');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])  
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['L5: r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end  
     set(gca,'Ydir','reverse');
xlabel('Fraction total input'); set(gca,'FontSize',10);xlim([0. 0.4]);%text(0.05,120,'L5')
end