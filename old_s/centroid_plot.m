function centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,gr,fe,label_name)

set(0,'DefaultLegendAutoUpdate','off')
if gr==0
    
ang2=out_ang_exL23;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 300]);
subplot(1,2,1);
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'m.','MarkerSize',10);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor='m';h4.MarkerEdgeColor='m';h4.MarkerSize=5;grid on;%h4.MarkerFaceAlpha=0.5
hold on;%text(-250,50,'L2/3');%legend('L2/3');legend boxoff
grid on;%h4.MarkerFaceAlpha=0.5
set(gca,'FontSize',12);
ang2=out_ang_exL4
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'g.','MarkerSize',10);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor='g';h4.MarkerEdgeColor='g';h4.MarkerSize=5;grid on;
set(gca,'Ydir','reverse');set(gca,'FontSize',12);
hold on;hold on;%text(-250,200,'L4');%legend(ang2,'L4');legend boxoff
ang2=out_ang_exL5
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'.','Color',[0.5 0.5 0.5],'MarkerSize',10);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor=[0.5 0.5 0.5];h4.MarkerEdgeColor=[0.5 0.5 0.5];h4.MarkerSize=5;grid on;
set(gca,'Ydir','reverse');hold on;line([0 0], [-4*69 8*69],'Color','k','LineStyle','-');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','-');hold on;plot(0,0,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on;title('EX','Color','r');ylabel('Vertical distance (µm)');xlabel('Horizontal distance (µm)');hold on;%text(-250,350,'L5');%legend('L5');legend boxoff
yticks([-200:200:600]);set(gca,'FontSize',12);
box off;

ang2=out_ang_inL23;
subplot(1,2,2);
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'m.','MarkerSize',10);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor='m';h4.MarkerEdgeColor='m';grid on;%h4.MarkerFaceAlpha=0.5
hold on;set(gca,'FontSize',10);%legend('L2/3');legend boxoff
ang2=out_ang_inL4;
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'g.','MarkerSize',9);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor='g';h4.MarkerEdgeColor='g';grid on;
set(gca,'Ydir','reverse');%legend('L4');legend boxoff
hold on;set(gca,'FontSize',12);
ang2=out_ang_inL5
h4 = plot(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,'.','Color',[0.5 0.5 0.5],'MarkerSize',10);xlim([-4*69 4*69]);ylim([-4*69 8*69]);%h4.MarkerFaceColor=[0.5 0.5 0.5];h4.MarkerEdgeColor=[0.5 0.5 0.5];grid on;
set(gca,'Ydir','reverse');hold on;line([0 0], [-4*69 8*69],'Color','k','LineStyle','-');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','-');hold on;plot(0,0,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on;title('IN','Color','b');xlabel('Horizontal distance (µm)');
box off;%legend('L5');legend boxoff
yticks([-200:200:600]);set(gca,'FontSize',12);
hold on;
plot(200,-200,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);text(220,-200,'L2/3');
hold on;
plot(200,-150,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);text(220,-150,'L4');
hold on;
plot(200,-100,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);text(220,-100,'L5');
grid on;
set(gca,'FontSize',12);

elseif gr==1
    [cmap]=inferno;
ang2=out_ang_exL23;    
pointsize=10;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 300]);
subplot(1,2,1);
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;

colormap(cmap);%text(-250,50,'L2/3');
ang2=out_ang_exL4;
yticks([-200:200:600])
hold on;
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled','s');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;%text(-250,150,'L4');

colormap(cmap);
ang2=out_ang_exL5
hold on;
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled','<');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;
set(gca,'Ydir','reverse');hold on;line([0 0], [-4*69 8*69],'Color','k','LineStyle','-');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','-');hold on;hold on;plot(0,0,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on;;title('EX Centre of mass');
c=colorbar;c.Label.String=label_name{1};c.Location='southoutside';  
ylabel('Vertical distance CoM (µm)');xlabel('Horizontal distance CoM (µm)')
%text(-250,350,'L5')
ang2=out_ang_inL23;
subplot(1,2,2);
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;

colormap(cmap);
ang2=out_ang_inL4;
hold on;
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled','s');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;

colormap(cmap);
ang2=out_ang_inL5
hold on;
h4 = scatter(ang2(a,3)*69-ang2(a,1)*69,ang2(a,4)*69-ang2(a,2)*69,pointsize,fe,'filled','<');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;
set(gca,'Ydir','reverse');hold on;line([0 0], [-4*69 8*69],'Color','k','LineStyle','-');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','-');hold on;hold on;plot(0,0,'^','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on;title('IN Centre of mass');c=colorbar;c.Label.String=label_name{1};c.Location='southoutside';    
yticks([-200:200:600]);
ylabel('Vertical distance CoM (µm)');xlabel('Horizontal distance CoM (µm)')
elseif gr==2
 
   pointsize=30;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);
 xin=out_ang_inL23(a,3)-out_ang_inL23(a,1); 
   xex=out_ang_exL23(a,3)-out_ang_exL23(a,1); 
   yin=out_ang_inL23(a,4)-out_ang_inL23(a,2); 
   yex=out_ang_exL23(a,4)-out_ang_exL23(a,2); 
    q=quiver(xex*69,yex*69,(xin-xex)*69,(yin-yex)*69,0);
   q.Color=[0.5 0.5 0.5];q.MaxHeadSize=0;q.LineWidth=1;
  hold on; h4 = scatter(xex*69,yex*69,pointsize,fe,'filled');xlim([-4*69 4*69]);ylim([-4*69 8*69]);grid on;
   set(gca,'Ydir','reverse');
   line([0 0], [-4*69 8*69],'Color','k','LineStyle','-');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','-');hold on;plot(0,0,'^r');
hold on;title('Centre of mass');c=colorbar;c.Label.String=label_name{1};    
[cmap]=buildcmap('ybk');colormap(cmap);box off
hold on;
  xex=out_ang_exL4(a,3)-out_ang_exL4(a,1); 
   yin=out_ang_inL4(a,4)-out_ang_inL4(a,2); 
   yex=out_ang_exL4(a,4)-out_ang_exL4(a,2); 
    q=quiver(xex*69,yex*69,(xin-xex)*69,(yin-yex)*69,0);
   q.Color=[0.5 0.5 0.5];q.MaxHeadSize=0;q.LineWidth=1;
 hold on; h4 = scatter(xex*69,yex*69,pointsize,fe,'filled')
 hold on;
  xex=out_ang_exL5(a,3)-out_ang_exL5(a,1); 
   yin=out_ang_inL5(a,4)-out_ang_inL5(a,2); 
   yex=out_ang_exL5(a,4)-out_ang_exL5(a,2); 
    q=quiver(xex*69,yex*69,(xin-xex)*69,(yin-yex)*69,0);
   q.Color=[0.5 0.5 0.5];q.MaxHeadSize=0;q.LineWidth=1;
 hold on; h4 = scatter(xex*69,yex*69,pointsize,fe,'filled')
xlabel('Horizontal distance (µm)');ylabel('Vertical distance (µm)')

else
    pointsize=30;
   fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);
 xin=out_ang_inL23(a,3)-out_ang_inL23(a,1); 
   xex=out_ang_exL23(a,3)-out_ang_exL23(a,1); 
   yin=out_ang_inL23(a,4)-out_ang_inL23(a,2); 
   yex=out_ang_exL23(a,4)-out_ang_exL23(a,2); 
    q=quiver(xex-xex,yex-yex,(xin-xex)*69,(yin-yex)*69,0);
   q.Color=[0.5 0.5 0.5];q.MaxHeadSize=0;q.LineWidth=1;
 % hold on; h4 = scatter(xin*69-xex,yin*69-yex,pointsize,fe,'filled');grid on;
   %set(gca,'Ydir','reverse'); 
   dx=(xin-xex)*69;
   dy=(yin-yex)*69;
   hypo=sqrt(((dx.^2)+(dy.^2)))
   sina=dy./hypo;
ang_a=asind(sina);
end


end