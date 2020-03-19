function corr_plot(par1,par2,par3,label_names)
if isempty(par3)==1
figure;set(gcf,'color','w');
 [R P]=corrcoef(par1,par2);scatter(par1,par2,25,'o','MarkerEdgeColor','k');box off;xlabel(label_names{1});ylabel(label_names{2});
%set(gca,'Ydir','reverse');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'g-');set(gca,'box','off');set(gcf,'color','w');axis square;   



else
[cmap]=buildcmap('kmgy');
pointsize=20;
figure; set(gcf,'color','w');
scatter(par1,par2,pointsize,par3,'filled');
box off; hold on; c=colorbar;colormap(cmap);c.Label.String = label_names{3}; %c.Ticks=[min(par3):round(min(par3)/3,-1):max(par3)];  
caxis([min(par3) max(par3)]);

if any(isnan(par2))==1
idxn = find(~isnan(par2));
else
idxn=1:length(par2);
end
P = polyfit(par1(idxn),par2(idxn),1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
  hold on;
  p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel(label_names{2});xlabel(label_names{1});
 [R P]=corrcoef(par1,par2,'row','complete');
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
end
end