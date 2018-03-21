%% Clean up
clearvars
close all
%% Load the files

%loading path (might have to edit this, depending on your own path)
load_path = 'R:\Share\Simon\Drago_Volker_Simon\Morphology\';

%file names (make sure the variable inside the file is named the same, I
%think "output" was different and I had to edit it)
file_names = {'morph_all'};

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
    all_cell{files} = all_cell{files}.morph_all;
    
end
%condense the files in a single cell
all_cell = vertcat(all_cell{:});

%get the cell number
cell_num = size(all_cell.totbrp,2);

%define the fields of interest
field_vec = [3:10,12:14,32:35,44];

%assemble a matrix with the relevant parameters
morpho_norm = zeros(cell_num,length(field_vec));
%get the structure field names
names = fieldnames(all_cell);

%for all the relevant fields
for fields = 1:length(field_vec)
    %load the matrix
    morpho_norm(:,fields) = all_cell.(names{field_vec(fields)});
end

%create a cell with the variable names
morpho_vars = names(field_vec);

%load the morpho cell names
morpho_names = all_cell.names;

%exclude the names of cells with less than 10 branching points
morpho_names = morpho_names(morpho_norm(:,1)>10);
cell_num = sum(morpho_norm(:,1)>10);
%standardize them
%for all the cells
for cells = 1:cell_num
    morpho_names{cells} = strcat(morpho_names{cells}(7:end-8),...
        morpho_names{cells}(end-3:end));
end

%exclude cells with less than 10 total branching points, since they are
%probably not filled enough
morpho_norm = morpho_norm(morpho_norm(:,1)>10,:);
%correct the pial distance
morpho_norm(:,16) = 1-morpho_norm(:,16);


%also save a copy in case I erase it to cluster something else
morpho_raw = morpho_norm;
%normalize the data
morpho_norm = zscore(morpho_norm);
%% Plot distributions of the parameters

close all

%get the number of parameters
paramt_num = size(morpho_norm,2);

figure
%for all the parameters
for paramt = 1:paramt_num
    subplot(round(sqrt(paramt_num)),ceil(sqrt(paramt_num)),paramt)
    histogram(morpho_raw(:,paramt))
    title(morpho_vars{paramt},'Interpreter','none')
end
%% Implement PCA

close all

%run a PCA on the normalized data
[morpho_coeff,morpho_score,morpho_latent] = pca(morpho_norm);

figure
imagesc(morpho_score)
ylabel('Cells')
xlabel('PC space variables')
title('PC-weighted data')

figure
plot(morpho_latent./sum(morpho_latent))
hold('on')
plot(cumsum(morpho_latent./sum(morpho_latent)))
title('PC normalized variance')

figure
imagesc(morpho_coeff)
set(gca,'YTick',1:size(morpho_norm,2),'YTickLabels',morpho_vars)
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
title('PCA loadings per variable')

figure
plot3(morpho_score(:,1),morpho_score(:,2),morpho_score(:,3),'*')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
%% Cluster the data

close all hidden

%define the number of morpho_clusters
morpho_clunum = 5;

%normalize the data and exclude cells with NaNs
% ep_cluster = morpho_norm;
% ep_cluster = ep_cluster(~isnan(sum(ep_cluster,2)),:);
% ep_cluster = normr_2(ep_cluster,2);
ep_cluster = morpho_score(:,1:2);

%cluster using linkage and cluster so that I can get the indexes
l_tree = linkage(ep_cluster,'average','euclidean');
morpho_clusters = cluster(l_tree,'maxclust',morpho_clunum);
%code to order the leaves more optimally
D = pdist(ep_cluster);
leafOrder = optimalleaforder(l_tree,D);
%plot the dendrogram
figure
dendrogram(l_tree,0,'Reorder',leafOrder,'Orientation','top')

%allocate memory to store the cluster averages in PC space
clu_ave = zeros(morpho_clunum,size(ep_cluster,2));
%and also for the averages and std in real space
clu_real = zeros(morpho_clunum,size(morpho_norm,2),2);
%save the number of members too
clu_mem = zeros(morpho_clunum,1);

%plot maps of a set of morpho_clusters
%for all the morpho_clusters
for clu = 1:morpho_clunum
    %average the cells in question
    clu_ave(clu,:) = squeeze(mean(ep_cluster(morpho_clusters==clu,:),1));
    %and also calc de std
    clu_real(clu,:,1) = squeeze(mean(morpho_norm(morpho_clusters==clu,:)));
    clu_real(clu,:,2) = squeeze(std(morpho_norm(morpho_clusters==clu,:))./sqrt(sum(morpho_clusters==clu)));
    %and count the number of members
    clu_mem(clu) = sum(morpho_clusters==clu);
end

figure
imagesc(clu_ave)
set(gca,'YTick',1:morpho_clunum,'YTickLabels',{clu_mem(:)})
set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',{'PC1','PC2','PC3'},'XTickLabelRotation',90)
set(gca,'TickLabelInterpreter','none')
ylabel('Number of cells')
colorbar
%also plot the cells sorted and label the different origins
figure
[~,ind] = sort(morpho_clusters);
% subplot(1,8,1:7)
imagesc(morpho_norm(ind,:));
ylabel('Cells')
% subplot(1,8,8)
% imagesc(ori_vec(ind))
% colormap(gca,[0 0 0;1 0 0;0 0 1])
%% Plot variables from the morpho_clusters

%define the target variables
tar_var = [1 2 3 11 16];

%get the number of plots
plot_num = length(tar_var);

figure
%for all the plots
for plots = 1:plot_num
    subplot(round(sqrt(plot_num)),ceil(sqrt(plot_num)),plots)
    errorbar(1:morpho_clunum,clu_real(:,tar_var(plots),1),clu_real(:,tar_var(plots),2),'*')
    set(gca,'XLim',[0 morpho_clunum + 1])
    xlabel(morpho_vars{tar_var(plots)})
end
%% Identify cells in common with the overlap data set

close all

%load the overlap data set 
%loading path (might have to edit this, depending on your own path)
load_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';
%load the most recent file from the overlap ones
listing = dir(strcat(load_path,'\*_invitro.mat'));
dates = datetime({listing.date});
[~,ind] = max(dates);
load_name = listing(ind).name;
%define the full loading path
load_full = fullfile(load_path,load_name);

%load the names of the import cells (special variable since it needed to be
%remapped)
import_names = load('R:\Share\Simon\Drago_Volker_Simon\names_input.mat','names_input');
import_names = import_names.names_input;
%remove the extra apostrophe at the end
%for all the cells
for cells = 1:length(import_names)
    import_names{cells} = strcat(import_names{cells}(1:end-1));
end

%cross reference the names
[~,ia,ib] = intersect(morpho_names,import_names);

%load the data from the import cells
import_mat = load(load_full,'invitro_norm');
import_mat = import_mat.invitro_norm;

%create averages of the common cell import properties using the local
%morpho_clusters
%get the number of common cells
com_num = length(ia);
%get a vector with the cluster indexes of the common cells
com_clusters = morpho_clusters(ia);
%and get the data for the corresponding cells
com_cells = import_mat(ib,:);

%allocate memory for the average profiles
com_ave = zeros(morpho_clunum,size(import_mat,2));
%for all the morpho_clusters
for clu = 1:morpho_clunum
    %skip the morpho_clusters that don't have representatives in the imported
    %cells
    if sum(com_clusters==clu)==0
        continue
    end
    %get the common cells for this cluster
    com_ave(clu,:) = mean(com_cells(com_clusters==clu,:));
end

%plot the averages
figure
imagesc(com_ave)
%% Calculate the cluster averages
close all
%allocate memory for the cluster average
morpho_cluave = zeros(morpho_clunum,size(morpho_norm,2));
%and for the number of cluster members
morpho_clumem = zeros(morpho_clunum,1);
%for all the clusters
for clu = 1:morpho_clunum
    %calculate the average of this cluster
    morpho_cluave(clu,:) = mean(morpho_norm(morpho_clusters==clu,:),1);
    %and the number of members
    morpho_clumem(clu) = sum(morpho_clusters==clu);
end

figure
imagesc(morpho_cluave)
%% Save the invivo calculations

%define the save path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_morpho.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(save_path,save_name),'morpho_raw','morpho_norm','morpho_score','morpho_clunum',...
    'morpho_clusters','morpho_coeff','morpho_latent','morpho_names',...
    'morpho_vars','morpho_cluave','morpho_clumem')
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
% %and the associated number of morpho_clusters
% morpho_clunum = clu_vec(bic_min);
% %get the indexes from the best model
% morpho_clusters = clu_cell{bic_min,2};
% 
% %allocate memory to store the cluster averages
% clu_ave = zeros(morpho_clunum,size(ep_cluster,2));
% %save the number of members too
% clu_mem = zeros(morpho_clunum,1);
% 
% %plot maps of a set of morpho_clusters
% %for all the morpho_clusters
% for clu = 1:morpho_clunum
% %     %if the cluster only contains 1 cell
% %     if sum(morpho_clusters==clu)==1
% %         %skip it
% %         continue
% %     end
%     %average the cells in question
%     clu_ave(clu,:) = squeeze(mean(ep_cluster(morpho_clusters==clu,:),1));
%     %and count the number of members
%     clu_mem(clu) = sum(morpho_clusters==clu);
% end
% 
% figure
% imagesc(clu_ave)
% set(gca,'YTick',1:morpho_clunum,'YTickLabels',{clu_mem(:)})
% set(gca,'XTick',1:size(ep_cluster,2),'XTickLabels',morpho_vars,'XTickLabelRotation',90)
% set(gca,'TickLabelInterpreter','none')
% colorbar
%% Calculate correlations between variables

close all

%calculate correlation matrix between variables (not cells)
[rho,pval] = corr(morpho_norm);

%show the correlation matrix
figure
imagesc(rho)
set(gca,'XTick',1:size(morpho_norm,2),'XTickLabels',morpho_vars,'XTickLabelRotation',90)
set(gca,'YTick',1:size(morpho_norm,2),'YTickLabels',morpho_vars)

set(gca,'TickLabelInterpreter','none')
title('Correlation between variables')
colorbar