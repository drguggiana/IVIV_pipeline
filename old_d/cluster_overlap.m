function cluster_overlap(clusters_1,clusters_2,link_type,dist_type)
%% Calculate overlap between types of clusters

%1)clusters_1 and clusters_2 are the indexes for the clusters, ordered so
%they match for the different cells
%2)link type is the type of linkage, normally use complete
%3)dist type is the type of distance metric, normally euclidean


%assemble a common vector with the cluster assignments from whole map and
%index based clusterings
all_cluster = cat(2,clusters_1,clusters_2);
%cluster the vector lolololol
%cluster using linkage and cluster so that I can get the indexes
all_tree = linkage(all_cluster,link_type,dist_type);
% all_clusters = cluster(all_tree,'maxclust',map_clunum);
% invitro_clusters = cluster(l_tree,'cutoff',cutoff_dist);
%code to order the leaves more optimally
all_D = pdist(all_cluster);
all_leafOrder = optimalleaforder(all_tree,all_D);
%plot the dendrogram
figure
% dendrogram(all_tree,0,'orientation','top','Reorder',all_leafOrder)

dendroplot(all_tree,all_leafOrder,map_dist,all_cluster,{'Whole map','Fractions'})