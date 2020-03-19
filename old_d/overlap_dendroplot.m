function overlap_dendroplot(l_tree,varargin)

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
% label_vec = varargin{2};

%sort the parameter matrix according to rows
[param_mat,idx] = sortrows(param_mat(:,2:4));
%transpose the matrix
param_mat = param_mat';

%create new figure
figure

%define the size of the subplot matrix
rows = 20;
columns = 20;
% %define the target coordinates of the dendrogram
% tar_rows = 1:5;
% tar_cols = 1:19;
% subplot_range(rows,columns,tar_rows,tar_cols);
% %plot the dendrogram
% dendrogram(l_tree,length(unique(param_mat(1,:))),'orientation','top','Reorder',idx);
% %label
% ylabel('Relative distance (a.u.)')
% %clear the ticks  on the x axis
% set(gca,'XTick',[],'TickLength',[0 0])


% %reorder the parameter matrix based on the cluster positions in the
% %dendrogram
% param_mat = param_mat(outperm,:)';
%% Plot the input clusters
%define the colormap
cmap = 'parula';
%define the target coordinates of the param plot
tar_rows = 6:10;
tar_cols = 1:19;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot the parameter matrix
imagesc(param_mat(1,:))
%label the plot
set(gca,'YTick',[],'TickLength',[0 0])
set(gca,'XTick',[],'XLim',[0 size(param_mat,2)+1])
colormap(gca,cmap)
ylabel('Input map clusters')

%plot the colorbar
%determine the number of clusters
clu_num = length(unique(param_mat(1,:)));
%define the target coordinates of the param plot
tar_rows = 6:10;
tar_cols = 20;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot it and label
imagesc((clu_num:-1:1)')
set(gca,'TickLength',[0 0])
set(gca,'XTick',[])
set(gca,'YTick',1:clu_num,'YTickLabels',clu_num:-1:1,'YAxisLocation','right')
% ylabel('A.U.')
colormap(gca,cmap)
%% Plot the morpho clusters
%define the colormap
cmap = 'jet';
%define the target coordinates of the param plot
tar_rows = 11:15;
tar_cols = 1:19;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot the parameter matrix (blanking the NaNs)
nan_mat = ~isnan(param_mat(2,:));
imagesc(param_mat(2,:),'AlphaData',nan_mat)
%label the plot
set(gca,'YTick',[],'TickLength',[0 0])
set(gca,'XTick',[],'XLim',[0 size(param_mat,2)+1])
colormap(gca,cmap)
ylabel('Morphological clusters')

%plot the colorbar
%determine the number of clusters
clu_num = length(unique(param_mat(2,~isnan(param_mat(2,:)))));
%define the target coordinates of the param plot
tar_rows = 11:15;
tar_cols = 20;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot it and label
imagesc((clu_num:-1:1)')
set(gca,'TickLength',[0 0])
set(gca,'XTick',[])
set(gca,'YTick',1:clu_num,'YTickLabels',clu_num:-1:1,'YAxisLocation','right')
% ylabel('A.U.')
colormap(gca,cmap)
%% Plot the invivo clusters
%define the colormap
cmap = 'jet';
%define the target coordinates of the param plot
tar_rows = 16:20;
tar_cols = 1:19;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot the parameter matrix (blanking the NaNs)
nan_mat = ~isnan(param_mat(3,:));
imagesc(param_mat(3,:),'AlphaData',nan_mat)
%label the plot
set(gca,'YTick',[],'TickLength',[0 0])
set(gca,'XTick',[],'XLim',[0 size(param_mat,2)+1])
colormap(gca,cmap)
ylabel('In vivo clusters')


%plot the colorbar
%determine the number of clusters
clu_num = length(unique(param_mat(3,~isnan(param_mat(3,:)))));
%define the target coordinates of the param plot
tar_rows = 16:20;
tar_cols = 20;
subplot_range(rows,columns,tar_rows,tar_cols);
%plot it and label
imagesc((clu_num:-1:1)')
set(gca,'TickLength',[0 0])
set(gca,'XTick',[])
set(gca,'YTick',1:clu_num,'YTickLabels',clu_num:-1:1,'YAxisLocation','right')
% ylabel('A.U.')
colormap(gca,cmap)

