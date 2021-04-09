
function sort_fractions_input(Lfr,vh,des)

if vh==1
[sv idxsort] = sort(Lfr(:,1),'ascend');
[cmap]=buildcmap('wmk');

% fig5=figure;set(gcf,'color','w');set(fig5, 'Position', [200, 0, 600,400]);
% imagesc(Lfr(idxsort,:))
% %pcolor(Lfr(idxsort,:))
% tmp=(Lfr(idxsort,1)-Lfr(idxsort,2))+1.5;
% colormap(cmap);c.Label.String = 'Pial depth (µm)';c=colorbar;
% hold on
% plot(tmp,1:147)

fig6=figure;set(gcf,'color','w');set(fig6, 'Position', [200, 0, 200, 225]);
x=[1:length(Lfr); 1:length(Lfr)]';
y=[Lfr(idxsort,2)*-1 Lfr(idxsort,3)]
bh=barh(x,y);
  bh(2).FaceColor=[0 0 1];
  bh(1).FaceColor=[1 0 0];
 set(gca,'Ydir','reverse')
%p1=plot(y(:,1),x(:,1),'r-o','LineWidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2);
%p1.Color(4) = 0.2;
%set(gca,'Ydir','reverse')
box off;
% hold on;
% p1=plot(y(:,2),x(:,2),'b-o','LineWidth',0.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2);
% p1.Color(4) = 0.2;
hold on;
tmp=(Lfr(idxsort,2)-Lfr(idxsort,3))
bh=barh(x(:,1),tmp);
ylabel('Cells sorted by pial depth');
xlabel('Fraction input');ylim([1 length(y)]);
yticks([0:20:147])
xlim([-1 1])
text(-0.5,-4,'EX','Color','r');
text(0.5,-4,'IN','Color','b');
title(des)
set(gca,'FontSize',10)
else
    
    
[sv idxsort] = sort(Lfr(:,1),'descend');
[cmap]=buildcmap('wmk');
fig6=figure;set(gcf,'color','w');set(fig6, 'Position', [200, 0, 500, 300]);
x=[1:147; 1:147]';
y=[Lfr(idxsort,1)*-1 Lfr(idxsort,2)]
bh=bar(x,y);
 bh(2).FaceColor=[0 0 1];
 bh(1).FaceColor=[1 0 0];
 box off;
hold on;
tmp=(Lfr(idxsort,1)-Lfr(idxsort,2))
plot(1:147,tmp*-1,'Color','k','LineWidth',1);
xlabel('Cells');
ylabel('Fraction input');ylim([1 length(y)]);
xticks([0:20:147])
ylim([-16 16])
title(des)
end
end