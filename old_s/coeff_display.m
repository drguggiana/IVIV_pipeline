%Display variance explained 
function coeff_display(coeff1,coeff2,bin_num,hbin_num)

if isempty(coeff2)==1
 
fig2 = figure;
set(fig2, 'Name', 'PCs coeff aligned maps');
set(fig2, 'Position', [0, 0, 800, 150]);set(gcf,'color','w')
set(gcf,'color','w')
[cmap]=buildcmap('kwg');
for i = 1:3
    subplot(1,3,i)
    imagesc(reshape(coeff1(:,i),16+bin_num,16+hbin_num));
    axis square;colorbar;colormap(cmap);caxis([-0.3 0.3]);
    set(gca,'xtick',[]);
set(gca,'ytick',[]);
end
c=colorbar;colormap(cmap);
c.Label.String = 'Coefficient';



    
else
fig2 = figure;
set(fig2, 'Name', 'PCs coeff ex in aligned maps');
set(fig2, 'Position', [0, 0, 800, 300]);set(gcf,'color','w')
set(gcf,'color','w')
[cmap]=buildcmap('kwg');
%[cmap]=buildcmap('kwb');
for i = 1:3
    subplot(2,3,i)
    imagesc(reshape(coeff1(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.3 0.3]);
    set(gca,'xtick',[]);
set(gca,'ytick',[]);
end
c=colorbar;colormap(cmap);
c.Label.String = 'Coefficient';
hold on;

for i = 1:3
    subplot(2,3,i+3)
    imagesc(reshape(coeff2(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.3 0.3]);
    set(gca,'xtick',[]);
set(gca,'ytick',[]);
end
c=colorbar;colormap(cmap);
c.Label.String = 'Coefficient';
    
end
end