%% Clean up
clearvars
close all
%% Define constants

%define the variable names
% var_names = {'VRest','Rin','Tau','Rheo','AP_f','AP_th','AP_amp','AP_hw','AP_slope','AHP','AP_adr','PialD'};
var_names = {'MinVm','PeakVm','InitVm','MaxVmSlope','HalfVm','Amp','MaxAHP',...
    'RiseT','FallT','BaseW','HalfW','FirstS','Vrest','tau','Rin','Sag','Rheo','MaxSF','PialD'};
%% Load the files

%loading path (might have to edit this, depending on your own path)
load_path = 'R:\Share\Simon\Drago_Volker_Simon\Ephys\';

%file names (make sure the variable inside the file is named the same, I
%think "output" was different and I had to edit it)
% file_names = {'Wildtype_all','WFS1positive_all','WFS1negative_all'};
file_names = {'cleaned_ephys'};

%get the number of files
file_num = length(file_names);
%allocate memory for the files
all_cell = cell(file_num,1);
%create a vector to keep track of the origin file
ori_cell = cell(file_num,1);
%for all the files
for files = 1:file_num
    %load the file and extract from the cell
    all_cell{files} = load(fullfile(load_path,file_names{files}));
    all_cell{files} = all_cell{files}.data;
%     %for the WFS files, eliminate columnes 13-15, since they don't involve
%     %this analysis
%     if files > 1
%         all_cell{files} = all_cell{files}(:,[1:12,16]);
%     end
%     %for all files, remove the Sag (fourth parameter)
%     all_cell{files} = all_cell{files}(:,[1:3,5:end]);
    %store the number of cells
    ori_cell{files} = zeros(size(all_cell{files},1),1)+files;
end
%concatenate all the data
all_data = cat(1,all_cell{:});
ori_vec = cat(1,ori_cell{:});
% %fix scale for pial distance 120-375 um (for this data set)
% all_data(:,13) = all_data(:,13);

%also save a copy in case I erase it to cluster something else
all_copy = all_data;

% %remove NaNs
% all_data(isnan(all_data)) = 0;
% %remove the cells that have NaNs
% nan_vec = isnan(sum(all_data,2));
% all_data = all_data(~nan_vec,:);
all_dataRaw = all_data;

pialD = all_data(:,19);
all_data = all_data(:,1:18);
%normalize the data
all_data = zscore(all_data);
%% Plot distributions of the parameters

close all

%get the number of parameters
paramt_num = size(all_data,2);

figure
%for all the parameters
for paramt = 1:paramt_num
    subplot(round(sqrt(paramt_num)),ceil(sqrt(paramt_num)),paramt)
    histogram(all_data(:,paramt))
    title(var_names{paramt},'Interpreter','none')
end
%% Implement PCA

close all

%run a PCA on the normalized data
[coeff,score,latent] = pca(all_data);

figure
imagesc(score)
ylabel('Cells')
xlabel('PC space variables')
title('PC-weighted data')

figure
plot(latent./sum(latent))
hold('on')
yyaxis right
plot(cumsum(latent)/sum(latent))
title('PC normalized variance')

figure
imagesc(coeff)
set(gca,'YTick',1:size(all_data,2),'YTickLabels',var_names)
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
title('PCA loadings per variable')

figure
plot3(score(:,1),score(:,2),score(:,3),'*')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
%% Cluster the data

close all hidden

%define the number of clusters
clu_num = 7;

%normalize the data and exclude cells with NaNs
% ep_cluster = all_data;
% ep_cluster = ep_cluster(~isnan(sum(ep_cluster,2)),:);
% ep_cluster = normr_2(ep_cluster,2);
ep_cluster = score(:,1:5);

%cluster using linkage and cluster so that I can get the indexes
%define the type of linkage
link_type = 'complete';
%define the distance to use
dist_type = 'euclidean';
l_tree = linkage(ep_cluster,link_type,dist_type);
% %use a monotone transformation to avoid inversions
% l_tree(:,3) = l_tree(:,3).^2;
clusters = cluster(l_tree,'maxclust',clu_num);
%code to order the leaves more optimally
D = pdist(ep_cluster);
leafOrder = optimalleaforder(l_tree,D);
%plot the dendrogram
figure
dendrogram(l_tree,0,'Reorder',leafOrder,'Orientation','top')

%allocate memory to store the cluster averages in PC space
clu_ave = zeros(clu_num,size(ep_cluster,2));
%the pial distance average
clu_pial = zeros(clu_num,2);
%and also for the averages and std in real space
clu_real = zeros(clu_num,size(all_data,2),2);
%save the number of members too
clu_mem = zeros(clu_num,1);

%plot maps of a set of clusters
%for all the clusters
for clu = 1:clu_num
    %average the cells in question
    clu_ave(clu,:) = squeeze(mean(ep_cluster(clusters==clu,:),1));
    %calculate the average and std pial distance
    clu_pial(clu,1) = mean(pialD(clusters==clu));
    clu_pial(clu,2) = std(pialD(clusters==clu))./sqrt(sum(clusters==clu));
    %and also calc de mean and std of the real variables
    clu_real(clu,:,1) = squeeze(mean(all_data(clusters==clu,:)));
    clu_real(clu,:,2) = squeeze(std(all_data(clusters==clu,:))./sqrt(sum(clusters==clu)));
    %and count the number of members
    clu_mem(clu) = sum(clusters==clu);
end

figure
imagesc(clu_ave)
set(gca,'YTick',1:clu_num,'YTickLabels',{clu_mem(:)})
set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',{'PC1','PC2','PC3'},'XTickLabelRotation',90)
set(gca,'TickLabelInterpreter','none')
ylabel('Number of cells')
colorbar
%also plot the cells sorted and label the different origins
figure
[~,ind] = sort(clusters);
subplot(1,8,1:7)
imagesc(all_data(ind,:));
ylabel('Cells')
subplot(1,8,8)
imagesc(ori_vec(ind))
colormap(gca,[0 0 0;1 0 0;0 0 1])
%% Evaluate the quality of the clusters

close all

%use evalclusters and the different criteria it provides
%define the methods
methods = {'CalinskiHarabasz','DaviesBouldin','gap','silhouette'};
%Define an anonymous function to calculate the clusters for the different
%criteria
h_cluster = @(Clu_mat,Clu_num)(clusterdata(Clu_mat,'linkage','complete','maxclust',Clu_num));
%get the number of methods
method_num = length(methods);
%define the vector numbers to compute
clunum_vec = 1:20;
figure
%for all the methods
for method = 4%1:method_num
    %run the method of choice
    switch method
        case 3
            eval_obj = evalclusters(ep_cluster,'linkage',methods{method},'Distance',dist_type,'KList',clunum_vec);
        case 4
            eval_obj = evalclusters(ep_cluster,h_cluster,methods{method},'Distance',dist_type,'KList',clunum_vec);
        case {1,2}
            eval_obj = evalclusters(ep_cluster,h_cluster,methods{method},'KList',clunum_vec);
    end
    subplot(round(sqrt(method_num)),ceil(sqrt(method_num)),method)
    plot(eval_obj)
    hold('on')
    
    clu_mems = zeros(7,1);
    for clu = 1:7
        clu_mems(clu) = sum(eval_obj.OptimalY==clu);
    end
end
%% Generate cluster shuffles and check clustering consistency

close all

%define the input matrix
input_mat = all_data;
%define the number of shuffles
shuff_num = 100;
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
    if N > 5
        %run a PCA on the normalized data
        [~,rand_score,~] = pca(random_cluster);
    end

    r_tree = linkage(rand_score(:,1:5),link_type,dist_type);
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
tree_vec = randperm(shuff_num,5);
figure
dendrogram(l_tree,0,'Reorder',leafOrder,'Orientation','top')
%get the y axis height for the original dendrogram
y_axis = get(gca,'YLim');
for trees = 1:5
    figure
    dendrogram(shuff_cell{trees,1},0,'Reorder',shuff_cell{trees,2},'Orientation','top')
    %apply a common y axis
    set(gca,'YLim',y_axis)
end

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
cutoff_prctile = prctile(concat_shuff,95);

%use the cutoff to cluster
clusters = cluster(l_tree,'cutoff',cutoff_prctile,'criterion','distance');
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
%% OFF Using the calculated cluster number, calculate the number of cluster members in a shuffle

% close all
% 
% %define the input matrix
% input_mat = all_data;
% %define the number of shuffles
% shuff_num = 1000;
% %allocate memory to store the results
% shuff_mems = zeros(shuff_num,clunum_prctile);
% %get the size of the clustering matrix
% [M,N] = size(input_mat);
% %for all the shuffles
% for shuffs = 1:shuff_num
%     
%     %randomize the clustering matrix along columns (across cells)
% %     random_cluster = ep_cluster();
%     
%     % Preserve the row indices
%     ColIndex = repmat((1:N),[M,1]);
%     % Get randomized column indices by sorting a second random array
%     [~,randomizedRowIndex] = sort(rand(M,N),1);
%     % Need to use linear indexing to create B
%     newLinearIndex = sub2ind([M,N],randomizedRowIndex,ColIndex);
%     random_cluster = input_mat(newLinearIndex);
%     
%     %if the matrix has more than the fice desired PCs, run PCA
%     if N > 5
%         %run a PCA on the normalized data
%         [~,rand_score,~] = pca(random_cluster);
%     end
% 
%     r_tree = linkage(rand_score(:,1:5),link_type,dist_type);
%     r_clusters = cluster(r_tree,'maxclust',clunum_prctile);
%     %code to order the leaves more optimally
% %     r_D = pdist(rand_score);
% %     r_leafOrder = optimalleaforder(r_tree,r_D);
%     %calculate the number of members per cluster
%     %for all the clusters
%     for clu = 1:clunum_prctile
%         shuff_mems(shuffs,clu) = sum(r_clusters==clu);
%     end
%  
% end
% 
% %sort the numbers
% shuff_sort = sort(shuff_mems,2);
% %plot the results
% figure
% histogram(shuff_sort(:,1))
% figure
% imagesc(shuff_mems)
% prctile(shuff_sort(:,1),95)
%% OFF Calculate the silhouette value for the calculated clusters
% %define the number of shuffles
% shuff_num = 100;
% %allocate memory to store the results
% shuff_sil = cell(shuff_num,1);
% 
% %for all the shuffles
% for shuffs = 1:shuff_num
%     % Preserve the row indices
%     ColIndex = repmat((1:N),[M,1]);
%     % Get randomized column indices by sorting a second random array
%     [~,randomizedRowIndex] = sort(rand(M,N),1);
%     % Need to use linear indexing to create B
%     newLinearIndex = sub2ind([M,N],randomizedRowIndex,ColIndex);
%     random_cluster = input_mat(newLinearIndex);
%     
%     %if the matrix has more than the fice desired PCs, run PCA
%     if N > 5
%         %run a PCA on the normalized data
%         [~,rand_score,~] = pca(random_cluster);
%     end
% 
%     r_tree = linkage(rand_score(:,1:5),link_type,dist_type);
%     r_clusters = cluster(r_tree,'maxclust',clunum_prctile);
%     h_cluster = @(Clu_mat,Clu_num)(clusterdata(Clu_mat,'linkage','complete','maxclust',Clu_num));
%     eval_obj = evalclusters(random_cluster,h_cluster,'silhouette','Distance',dist_type,'KList',clunum_prctile);
%     shuff_sil{shuffs} = eval_obj.ClusterSilhouettes{1};
% end
%% OFF Plot silhouette shuffle results

% close all
% 
% figure
% %get all the values combined
% all_sil = sort(horzcat(shuff_sil{:}),1);
% histogram(all_sil)
% % boxplot(all_sil')
% prctile(all_sil(all_sil~=1),95)
% 
% %show the original results
% eval_obj = evalclusters(ep_cluster,h_cluster,'silhouette','Distance',dist_type,'KList',clunum_prctile);
% 
% %for all the clusters
% %allocate memory to store the cluster averages in PC space
% clu_ave = zeros(clunum_prctile,size(ep_cluster,2));
% 
% %save the number of members too
% clu_label = cell(clunum_prctile,1);
% 
% %plot maps of a set of clusters
% %for all the clusters
% for clu = 1:clunum_prctile
%     %average the cells in question
%     clu_ave(clu,:) = squeeze(mean(ep_cluster(eval_obj.OptimalY==clu,:),1));
%     
%     %and count the number of members
%     clu_label{clu} = strcat(num2str(sum(eval_obj.OptimalY==clu)),'_',num2str(eval_obj.ClusterSilhouettes{1}(clu)));
% end
% 
% figure
% imagesc(clu_ave)
% set(gca,'YTick',1:clunum_prctile,'YTickLabels',clu_label)
% set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',{'PC1','PC2','PC3'},'XTickLabelRotation',90)
% set(gca,'TickLabelInterpreter','none')
% ylabel('Number of cells')
% colorbar
%% OFF Compute cluster members distance to cluster center in a shuffle 

% close all
% 
% %define the input matrix
% input_mat = all_data;
% %define the number of shuffles
% shuff_num = 100;
% %allocate memory to store the results
% shuff_cludist = zeros(shuff_num,clunum_prctile);
% %get the size of the clustering matrix
% [M,N] = size(input_mat);
% %for all the shuffles
% for shuffs = 1:shuff_num
%     
%     %randomize the clustering matrix along columns (across cells)
%     
%     % Preserve the row indices
%     ColIndex = repmat((1:N),[M,1]);
%     % Get randomized column indices by sorting a second random array
%     [~,randomizedRowIndex] = sort(rand(M,N),1);
%     % Need to use linear indexing to create B
%     newLinearIndex = sub2ind([M,N],randomizedRowIndex,ColIndex);
%     random_cluster = input_mat(newLinearIndex);
%     
%     %if the matrix has more than the fice desired PCs, run PCA
%     if N > 5
%         %run a PCA on the normalized data
%         [~,rand_score,~] = pca(random_cluster);
%     end
% 
%     r_tree = linkage(rand_score(:,1:5),link_type,dist_type);
%     r_clusters = cluster(r_tree,'maxclust',clunum_prctile);
%     %code to order the leaves more optimally
% %     r_D = pdist(rand_score);
% %     r_leafOrder = optimalleaforder(r_tree,r_D);    
%     
%     %save the results
%     %for all the clusters
%     for clu = 1:clunum_prctile
%         %get the number of members
%         shuff_nummems = sum(r_clusters==clu);
%         %if it's a singleton cluster, assign a 0 to it
%         if shuff_nummems == 1
%             continue
%         end
%         %get the members of the cluster
%         shuff_clumems = rand_score(r_clusters==clu,1:5);
%         
%         %get the cluster average
%         shuff_ave = mean(shuff_clumems);
%         %get the distances between the average and the members of the
%         %cluster
%         shuff_D = pdist([shuff_ave;shuff_clumems]);
%         %calculate the average of the distances to the average
%         shuff_cludist(shuffs,clu) = mean(shuff_D(1:shuff_nummems-1));
%     end
%     
% end
% 
% figure
% histogram(shuff_cludist(:))
% 
% prctile(shuff_cludist(:),95)
% 
% %calculate the distances in the real data
%% Plot variables from the clusters

%define the target variables
tar_var = [5 6 3 4 1 19];

%get the number of plots
plot_num = length(tar_var);

figure
%for all the plots
for plots = 1:plot_num
    subplot(round(sqrt(plot_num)),ceil(sqrt(plot_num)),plots)
    if tar_var(plots) == 19
        errorbar(1:clu_num,clu_pial(:,1),clu_pial(:,2),'*')
    else
        errorbar(1:clu_num,clu_real(:,tar_var(plots),1),clu_real(:,tar_var(plots),2),'*')
    end
    set(gca,'XLim',[0 clu_num + 1])
    xlabel(var_names{tar_var(plots)})
end
%% OFF Try GMM

% close all
% 
% %define a vector with cluster numbers
% clu_vec = [2 3 4 5 6 8 10];
% %get the number of cluster runs
% clu_size = length(clu_vec);
% %create a cell array to save the clustering results
% clu_cell = cell(clu_size,3);
% %define the statistical settings
% opts = statset('MaxIter',1000);
% %for all the cluster numbers
% for clu = 1:clu_size
%     
%     %create a GMM for the data
%     clu_cell{clu,1} = fitgmdist(ep_cluster,clu_vec(clu),'RegularizationValue',0.0001,...
%         'CovarianceType','diagonal','Replicates',100,'Options',opts);
%     
%     %cluster the data accordingly
%     clu_cell{clu,2} = cluster(clu_cell{clu,1},ep_cluster);
%     
%     %and get the BIC value
%     clu_cell{clu,3} = clu_cell{clu,1}.BIC;
% end
% 
% %plot the BIC
% figure
% plot(clu_vec,vertcat(clu_cell{:,3}));
% 
% %get the BIC minimum coordinate
% [~,bic_min] = min(vertcat(clu_cell{:,3}));
% %and the associated number of clusters
% clu_num = clu_vec(bic_min);
% %get the indexes from the best model
% clusters = clu_cell{bic_min,2};
% 
% %allocate memory to store the cluster averages
% clu_ave = zeros(clu_num,size(ep_cluster,2));
% %save the number of members too
% clu_mem = zeros(clu_num,1);
% 
% %plot maps of a set of clusters
% %for all the clusters
% for clu = 1:clu_num
% %     %if the cluster only contains 1 cell
% %     if sum(clusters==clu)==1
% %         %skip it
% %         continue
% %     end
%     %average the cells in question
%     clu_ave(clu,:) = squeeze(mean(ep_cluster(clusters==clu,:),1));
%     %and count the number of members
%     clu_mem(clu) = sum(clusters==clu);
% end
% 
% figure
% imagesc(clu_ave)
% set(gca,'YTick',1:clu_num,'YTickLabels',{clu_mem(:)})
% set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',var_names,'XTickLabelRotation',90)
% set(gca,'TickLabelInterpreter','none')
% colorbar
%% Calculate correlations between variables

close all

%calculate correlation matrix between variables (not cells)
[rho,pval] = corr(all_dataRaw);

%show the correlation matrix
figure
imagesc(rho)
set(gca,'XTick',1:size(all_dataRaw,2),'XTickLabels',var_names,'XTickLabelRotation',90)
set(gca,'YTick',1:size(all_dataRaw,2),'YTickLabels',var_names)

set(gca,'TickLabelInterpreter','none')
title('Correlation between variables')
colorbar