function [shuff_cell,clunum_prctile,cutoff_dist] = shuffle_cluster_1(input_mat,...
    shuff_num,pc_num,ep_cluster,link_type,dist_type,cutoff_prctile)

%% Generate cluster shuffles and check clustering consistency

close all

% %define the input matrix
% input_mat = all_data;
% %define the number of shuffles
% shuff_num = 100;
%allocate memory to store the results
shuff_cell = cell(shuff_num,3);
%get the size of the clustering matrix
[M,N] = size(input_mat);
%for all the shuffles
for shuffs = 1:shuff_num
    
    %randomize the clustering matrix along columns (across cells)
    
    % Preserve the row indices
    ColIndex = repmat((1:N),[M,1]);
    % Get randomized column indices by sorting a second random array
    [~,randomizedRowIndex] = sort(rand(M,N),1);
    % Need to use linear indexing to create B
    newLinearIndex = sub2ind([M,N],randomizedRowIndex,ColIndex);
    random_cluster = input_mat(newLinearIndex);
    
    %if the matrix has more than the fice desired PCs, run PCA
    if N > pc_num
        %run a PCA on the normalized data
        [~,rand_score,~] = pca(random_cluster);
    else
        rand_score = random_cluster;
    end

    r_tree = linkage(rand_score(:,1:pc_num),link_type,dist_type);
    %code to order the leaves more optimally
    r_D = pdist(rand_score);
    r_leafOrder = optimalleaforder(r_tree,r_D);    
    
    %save the results
    shuff_cell{shuffs,1} = r_tree;
    shuff_cell{shuffs,2} = r_leafOrder;
    shuff_cell{shuffs,3} = r_D;
    
end
%% Plot the shuffle results
close all
%pick 5 trees at random and display
% tree_vec = randperm(shuff_num,5);

%calculate the original distance tree
l_tree = linkage(ep_cluster,link_type,dist_type);
%code to order the leaves more optimally
D = pdist(ep_cluster);
leafOrder = optimalleaforder(l_tree,D);
figure
dendrogram(l_tree,0,'Reorder',leafOrder,'Orientation','top')
%get the y axis height for the original dendrogram
y_axis = get(gca,'YLim');
% for trees = 1:5
%     figure
%     dendrogram(shuff_cell{trees,1},0,'Reorder',shuff_cell{trees,2},'Orientation','top')
%     %apply a common y axis
%     set(gca,'YLim',y_axis)
% end

%Average the distance vectors and plot a tree
ave_D = mean(vertcat(shuff_cell{:,3}),1);
ave_tree = linkage(ave_D,link_type);
ave_leaforder = optimalleaforder(ave_tree,ave_D);
figure
dendrogram(ave_tree,0,'Reorder',ave_leaforder,'Orientation','top')
%apply a common y axis
set(gca,'YLim',y_axis)

%plot the distributions of the real distances and the shuffled ones
figure
concat_shuff = horzcat(shuff_cell{:,3});
histogram(concat_shuff,'Normalization','pdf')
hold('on')
histogram(D,'Normalization','pdf')
%get the 95th percentile for a cutoff
cutoff_dist = prctile(concat_shuff,cutoff_prctile);


%use the cutoff to cluster
clusters = cluster(l_tree,'cutoff',cutoff_dist,'criterion','distance');
%get the number of clusters
clunum_prctile = size(unique(clusters),1);
%allocate memory to store the cluster averages in PC space
clu_ave = zeros(clunum_prctile,size(ep_cluster,2));
% %the pial distance average
% clu_pial = zeros(clu_num,2);
% %and also for the averages and std in real space
% clu_real = zeros(clu_num,size(all_data,2),2);
%save the number of members too
clu_mem = zeros(clunum_prctile,1);

%plot maps of a set of clusters
%for all the clusters
for clu = 1:clunum_prctile
    %average the cells in question
    clu_ave(clu,:) = squeeze(mean(ep_cluster(clusters==clu,:),1));
    
%     %and also calc de mean and std of the real variables
%     clu_real(clu,:,1) = squeeze(mean(all_data(clusters==clu,:)));
%     clu_real(clu,:,2) = squeeze(std(all_data(clusters==clu,:))./sqrt(sum(clusters==clu)));
    %and count the number of members
    clu_mem(clu) = sum(clusters==clu);
end

figure
imagesc(clu_ave)
set(gca,'YTick',1:clunum_prctile,'YTickLabels',{clu_mem(:)})
set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',{'PC1','PC2','PC3'},'XTickLabelRotation',90)
set(gca,'TickLabelInterpreter','none')
ylabel('Number of cells')
colorbar