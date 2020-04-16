function dendroplot(l_tree,leafOrder,cutoff_dist,varargin)

%this function plots a dendrogram with the supplied hierarchical tree
%information, colors the clusters in it and then also displays an arbitrary
%matrix of parameters underneath, followed by the soma depths of the cells
%provided

%Parameter list
%1) Distance tree obtained from the linkage function
%2) Optimize leaf order for the hierarchical tree, obtained from the
%optimaleaforder function
%3) Cutoff distance to determine the tree clusters (as determined from the
%cluster shuffling script for example)
%4) Optional: parameter matrix to plot under the dendrogram
%(cells,variables)
%5) Optional: label vector for the variables in the parameter matrix
%6) Optional: matrix with the x and y coordinates of each cell's soma

%assign the variables based on the input
param_mat = varargin{1};
label_vec = varargin{2};
color_idx={'m','g','c','y'};

%[cmap]=buildcmap('kby');
[cmap]=buildcmap('wb')
%create new figure
figure
 set(gcf,'color','w');
%define the size of the subplot matrix
rows = 20;
columns = 20;
%define the target coordinates of the dendrogram
tar_rows = 1:10;
tar_cols = 1:19;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot the dendrogram
[H,~,outperm] = dendrogram(l_tree,0,'orientation','top','Reorder',leafOrder,'ColorThreshold',cutoff_dist);
%label
ylabel('Relative distance (a.u.)')
%clear the ticks  on the x axis
set(gca,'XTick',[],'TickLength',[0 0])
set(H,'LineWidth',1.5);
%  lineColours = cell2mat(get(H,'Color'));
% for m=1:clus_numb
% idx_c=find(lineColours(:,m)==1);
% hold on;
% for i=1:length(idx_c)
%     H(idx_c(i)).Color=color_idx{m};
% end
% hold on;
% idx_c=[];
% % idx2=find(lineColours(:,2)==1);
% % idx1=find(lineColours(:,3)==1);
% end
 
%if the parameter matrix is supplied
if nargin > 3
    %reorder the parameter matrix based on the cluster positions in the
    %dendrogram
    param_mat = param_mat(outperm,:)';

    %define the target coordinates of the param plot
    tar_rows = 11:16;
    tar_cols = 1:19;
    subplot_range(rows,columns,tar_rows,tar_cols);
    %plot the parameter matrix
   b= imagesc(param_mat)
   colormap(cmap)
  
 set(b,'AlphaData',~isnan(param_mat))
    %label the plot
    set(gca,'YTick',1:size(param_mat,1),'YTickLabels',label_vec,'TickLength',[0 0])
    set(gca,'XTick',[],'XLim',[0 size(param_mat,2)+1])
    
    %plot the colorbar
    %define the target coordinates of the param plot
    tar_rows = 11:16;
    tar_cols = 20;
    subplot_range(rows,columns,tar_rows,tar_cols);
    %plot it and label
    imagesc((255:-1:0)')
    colormap(cmap)
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',[])
    set(gca,'YTick',[1 256],'YTickLabels',[1 0],'YAxisLocation','right')
    ylabel('A.U.')
    %axis off;
    box off
end

%if the soma depth is supplied
if nargin > 5
    pialD = varargin{3};
    %also plot the pial distance
    tar_rows = 17:20;
    tar_cols = 1:19;
    subplot_range(rows,columns,tar_rows,tar_cols);
    %plot the soma depths and label
    %plot(1:size(param_mat,2),pialD(outperm,2),'k^','MarkerFaceColor','k')
     plot(1:size(param_mat,2),pialD(outperm),'k^','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',3)
     set(gca,'YLim',[140 360],'TickLength',[0 0],'YDir','reverse')
    set(gca,'XLim',[0 size(param_mat,2)+1],'XTick',[])
     ylabel('Soma Depth (um)')
     box off;axis off

end

