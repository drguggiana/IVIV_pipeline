%% Script to generate a figure of the cluster overlaps between data sets
%% Clean up
clearvars
close all
%% Load the data

%define the path
load_path = 'R:\Share\Simon\Drago_Volker_Simon\Cluster_str\clu_str.mat';

%load the structure
data = load(load_path);
data = data.clu_str;

%load the distance tree
l_tree = data.l_tree_input;
%% Create the matrices that include all cells

%allocate memory for the matrix
overlap_mat = zeros(size(data.input,1),length(fieldnames(data))+1);

%load the clusters for the input
overlap_mat(:,[1 2]) = data.input;

%for all the cells in the input map
for cells = 1:size(overlap_mat,1)
    %find the morpho cluster, otherwise NaN the position
    morpho_idx = data.morpho==overlap_mat(cells,1);
    
    if sum(morpho_idx)==0
        overlap_mat(cells,3) = NaN;
    else
        overlap_mat(cells,3) = data.morpho(morpho_idx,2);
    end
    
    %do the same with the invivo indexes
    invivo_idx = data.invivo==overlap_mat(cells,1);
    
    if sum(invivo_idx)==0
        overlap_mat(cells,4) = NaN;
    else
        overlap_mat(cells,4) = data.invivo(invivo_idx,2);
    end
end

% %exclude singlets
% %get the cluster ids for the input maps
% input_clu = unique(overlap_mat(:,2));
% %go through all the clusters
% for clu = input_clu'
%     %get the row ids corresponding to this cluster
%     clu_id = overlap_mat(:,2)==clu;
%     %if it is a singlet
%     if sum(clu_id)==1
%         %eliminate the entire row from the matrix
%         overlap_mat = overlap_mat(~clu_id,:);
%     end
% end

%flip the indexes of input clusters 1 and 3
clu_3 = overlap_mat(:,2)==3;
clu_1 = overlap_mat(:,2)==1;

overlap_mat(clu_1,2)= 3;
overlap_mat(clu_3,2)= 1;
%% Plot the overall figure

close all


link_type = 'complete';
dist_type = 'euclidean';
invitro_clunum = 3;
leaf_nodes = 3;
% l_tree = linkage(overlap_mat(:,2),link_type,dist_type);
% invitro_clusters = cluster(l_tree,'maxclust',invitro_clunum);
% invitro_clusters = cluster(l_tree,'cutoff',cutoff_dist);
%code to order the leaves more optimally
% D = pdist(overlap_mat(:,2));
% leafOrder = optimalleaforder(l_tree,D);

overlap_dendroplot(l_tree,overlap_mat)
% [~,~,outperm] = dendrogram(l_tree,leaf_nodes,'orientation','top','Reorder',leafOrder);
