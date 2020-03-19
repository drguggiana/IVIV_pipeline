
function display_sortfr(Lfr,vh,des)

if vh==1
[sv idxsort] = sort(Lfr(:,1),'descend');
[cmap]=buildcmap('wmk');

% fig5=figure;set(gcf,'color','w');set(fig5, 'Position', [200, 0, 600,400]);
% imagesc(Lfr(idxsort,:))
% %pcolor(Lfr(idxsort,:))
% tmp=(Lfr(idxsort,1)-Lfr(idxsort,2))+1.5;
% colormap(cmap);c.Label.String = 'Pial depth (µm)';c=colorbar;
% hold on
% plot(tmp,1:147)

fig6=figure;set(gcf,'color','w');set(fig6, 'Position', [200, 0, 300, 500]);
x=[1:147; 1:147]';
y=[Lfr(idxsort,1)*-1 Lfr(idxsort,2)]
bh=barh(x,y);
 bh(2).FaceColor=[0 0 1];
 bh(1).FaceColor=[1 0 0];
set(gca,'Ydir','reverse')
box off;
hold on;
tmp=(Lfr(idxsort,1)-Lfr(idxsort,2))
plot(tmp*-1,1:147,'Color','k','LineWidth',2);
ylabel('Cells');
xlabel('Fraction input');ylim([1 length(y)]);
yticks([0:20:147])
xlim([-0.4 0.4])
text(-0.2,-3,'EX','Color','r');)
text(0.2,-3,'IN','Color','b');
title(des)
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
plot(1:147,tmp*-1,'Color','k','LineWidth',2);
xlabel('Cells');
ylabel('Fraction input');ylim([1 length(y)]);
xticks([0:20:147])
ylim([-16 16])
title(des)
end
end