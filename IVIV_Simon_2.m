%% Clean up
clearvars
close all
%% Define constants

%define the layers
layer_c = {[1],[2 3 4 5],[6 7 8],9:16};
%get the number of layers
layer_num = length(layer_c);
%define the layer names
layer_name = {'L1','L2/3','L4','DL'};
%% Load the files

%loading path (might have to edit this, depending on your own path)
load_path = 'R:\Share\Simon\Drago_Volker_Simon\';

%file names (make sure the variable inside the file is named the same, I
%think "output" was different and I had to edit it)
file_names = {'all_invitro','invivo_all','output'};

%get the number of files
file_num = length(file_names);
%allocate memory for the files
all_cell = cell(file_num,1);
%for all the files
for files = 1:file_num
    %load the file and extract from the cell
    all_cell{files} = load(fullfile(load_path,file_names{files}));
    all_cell{files} = all_cell{files}.(file_names{files});
    
end
%get the names of the cells
cell_names = squeeze(all_cell{1}.names);
%and the corrected names, for saving(special variable since it needed to be
%remapped)
invitro_names = load('R:\Share\Simon\Drago_Volker_Simon\names_input.mat','names_input');
invitro_names = invitro_names.names_input;
%extract the extra apostrophe at the end of the file names
%for all the names
for names = 1:length(invitro_names)
    invitro_names{names} = invitro_names{names}(1:end-1);
end
%% Load and apply the cell flip for when the slide was turned over

%define the file path
flip_path = 'R:\Share\Simon\Drago_Volker_Simon\flip.mat';

%load the variable
load(flip_path)

%apply the flipping of the maps that should be flipped (1 means flip)
%for all the cells
for cells = 1:size(all_cell{1}.name,1)
    %skip the 0s
    if flip(cells) == 0
        continue
    end
    all_cell{1}.integ_epsc(:,:,cells) = all_cell{1}.integ_epsc(:,16:-1:1,cells);
    all_cell{1}.integ_ipsc(:,:,cells) = all_cell{1}.integ_ipsc(:,16:-1:1,cells);
end
%% Include the soma center data (mostly the vertical position)

close all
%get the soma center information
soma_cent = squeeze(all_cell{1}.soma)';

% %plot the distribution of positions
% figure
% histogram(soma_cent(:,1));
% title('X position distribution')
% 
% figure
% histogram(soma_cent(:,2));
% title('Y position distribution')

%also get the soma center of the cells with functional info
soma_func = squeeze(all_cell{3}.soma)';
%% Define the field of interest and collapse the e and i maps

%define the all cells field to work with
all_c = all_cell{1};

%load the maps into a single matrix
%get the number of vertical locations
v_num = size(all_c.integ_epsc,1);
%get the number of horizontal locations
h_num = size(all_c.integ_epsc,2);
%get the number of cells
cell_num = size(all_c.integ_epsc,3);
%allocate memory for the matrix
maps_mat = zeros(v_num,h_num,cell_num,2);

%collapse the maps into the matrix
maps_mat(:,:,:,1) = cat(3,all_c.integ_epsc);
maps_mat(:,:,:,2) = cat(3,all_c.integ_ipsc);

%calculate normalized versions
%allocate memory for the normalized version
norm_maps = zeros(v_num,h_num,cell_num,2);
for cells = 1:cell_num
    norm_maps(:,:,cells,1) = normr_2(maps_mat(:,:,cells,1)); 
    norm_maps(:,:,cells,2) = normr_2(maps_mat(:,:,cells,2)); 
end
%remove the NaNs (coming from maps that have only zeros)
norm_maps(isnan(norm_maps(:))) = 0;
%% OFF Generate overlap maps by binarizing the images
% close all
% %create the matrix with the maps
% bin_map = maps_mat>0;
% 
% %define the smooth factor
% sf = 10;
% %generate a blank matrix to fill in the other color channels
% blank = ones(size(bin_map,1),size(bin_map,2));
% blankw = ones(size(bin_map,1),size(bin_map,2),3);
% %calculate the average across cells for both E and I
% ave_map = squeeze(mean(bin_map,3));
% 
% % %set the excitation map as a red image, smoothing by sf
% exc_map = imresize(cat(3,blank,1-normr_2(ave_map(:,:,1)),1-normr_2(ave_map(:,:,1))),sf);
% % %and the inhibition map as a blue one
% inh_map = imresize(cat(3,1-normr_2(ave_map(:,:,2)),1-normr_2(ave_map(:,:,2)),blank),sf);
% % %blend the two images using alpha (and make it double cause default is 8bit
% im_ex = double(imfuse(exc_map,inh_map,'method','blend'));
% 
% figure
% 
% %normalize the image
% im_ex = normr_2(im_ex);
% %remove the NaNs (if there was an empty channel)
% im_ex(isnan(im_ex)) = 0;
% 
% %actually plot the map
% imagesc(im_ex)
%% OFF Also calculate overlap index based on the bin maps
% close all
% %allocate memory for the index
% ov_index = zeros(cell_num,1);
% 
% %for all the cells
% for cells = 1:cell_num
%     %get the overlap map
%     ov_map = squeeze(sum(bin_map(:,:,cells,:),4));
%     %quantify the index
%     ov_index(cells) = sum(ov_map(:)>1)./sum(ov_map(:)>0);
% end
% 
% %plot the distribution of indexes
% figure
% histogram(ov_index)
% title('Overlap Index')
%% Calculate overlap per layer

close all hidden

%create the matrix with the maps
bin_map = maps_mat>0;
%allocate memory for the indexes per layer
ov_perlayer = zeros(cell_num,layer_num);
%and the connection strength
% cs_perlayer = zeros(cell_num,layer_num);
cs_all = zeros(cell_num,1);

%for all the cells
for cells = 1:cell_num
    %get the overlap map
    ov_map = squeeze(sum(bin_map(:,:,cells,:),4));
    %and the diff map
%     diff_map = squeeze(diff(maps_mat(:,:,cells,:),1,4).^2);
    diff_map = sqrt(sum(maps_mat(:,:,cells,:).^2,4));
    %calculate the cell's total response "magnitude"
%     total_mag = sum(diff_map(:));
    cs_all(cells) = sum(diff_map(:));
    %for all the layers
    for layers = 1:layer_num

        %get the layers of interest
        ov_layers = ov_map(layer_c{layers},:);
        
        %get the index with the corresponding layers
        ov_perlayer(cells, layers) = sum(ov_layers(:)>1)./sum(ov_layers(:)>0);
%         %calculate the difference in mag for this layer and cell
%         layer_diff = diff_map(layer_c{layers},:);
        %scale the overlap index by the difference
%         cs_perlayer(cells, layers) = sum(layer_diff(:))/total_mag;
%         cs_perlayer(cells, layers) = total_mag;
        
    end
end

%remove NaNs
ov_perlayer(isnan(ov_perlayer)) = 0;

%plot a map of cells and layers
figure
imagesc(ov_perlayer')
set(gca,'YTick',1:4,'YTickLabels',layer_name,'XTick',[])
xlabel('Cells')
colorbar
title('Overlap per layer and cell')
%also plot the average ov per layer
figure
errorbar(1:4,nanmean(ov_perlayer,1),nanstd(ov_perlayer,0,1)./sqrt(cell_num),'*')
title('Overlap per layer')
set(gca,'XTick',1:4,'XTickLabels',layer_name,'XLim',[0 5])

%plot histograms per layer
figure
for layers = 1:layer_num
    subplot(2,2,layers)
    histogram(ov_perlayer(:,layers))
    title(layer_name{layers})
%     hold('on')
end
% %also plot the average cs per layer
% figure
% errorbar(1:4,nanmean(cs_perlayer,1),nanstd(cs_perlayer,0,1)./sqrt(cell_num),'*')
% title('Overlap per layer')
% set(gca,'XTick',1:4,'XTickLabels',layer_name,'XLim',[0 5])
% 
% %plot histograms of cs per layer
% figure
% for layers = 1:layer_num
%     subplot(2,2,layers)
%     histogram(cs_perlayer(:,layers))
%     title(layer_name{layers})
% %     hold('on')
% end
%% Load the functional info

%load it from the main cell
func_info = all_cell{2};
%change the responder property of cells that were misleadingly assigned as
%non-responder
mislead_vec = [4 6 16 23 31 35 50 54 67];
%replace the responder value
for cells = mislead_vec
    func_info(cells).responder = 1;
end
%select the parameter to plot with the invitro_clusters
%load the parameter of interest and make it binary
tar_par = vertcat(func_info.responder)==1;
% tar_par = vertcat(func_info.ave_responder)==1;
% tar_par = vertcat(func_info.ODscore_mean)>0;
% tar_par = vertcat(func_info.ODI_final)>0;
% tar_par = vertcat(func_info.OSI_contra)>0.5;
% tar_par = vertcat(func_info.OSI_ipsi)>0.5;
% tar_par = vertcat(func_info.DSI_contra)>0.5; %maybe something here
% tar_par = vertcat(func_info.DSI_ipsi)>0.5;
% tar_par = vertcat(func_info.responder_contra_n)==1;
% tar_par = vertcat(func_info.responder_ipsi_n)==1;

%non responder vector
% tar_par = zeros(size(tar_par));
% tar_par([7 8 12 20 22 24 27 29 45 46 47 60 62]) = 1;
% tar_par = tar_par==1;

% tar_par = zeros(143,1);
% tar_par(1:48) = 1;
% tar_par = tar_par==1;
%% OFF Normalize data for clustering
% close all hidden
% 
% %make NaNs 0s
% invitro_norm = normr_2([ov_perlayer,cs_all],2);
% invitro_norm(isnan(ov_perlayer)) = 0;
% 
% % ov_cluster = ov_cluster(1:48,:);
% 
% %add the normalized depth information
% soma_norm = normr_2(soma_cent(:,2));
% invitro_norm = cat(2,invitro_norm,soma_norm);
%% OFF PCA the data

% close all
% %define the variable names
% % invitro_vars = {'OL1','OL2/3','OL4','ODL','CL1','CL2/3','CL4','CDL','PD'};
% invitro_vars = {'OL1','OL2/3','OL4','ODL','CS','PD'};
% 
% %run a PCA on the normalized data
% [coeff,score,latent] = pca(invitro_norm);
% 
% figure
% imagesc(score)
% ylabel('Cells')
% xlabel('PC space variables')
% title('PC-weighted data')
% 
% figure
% plot(latent./sum(latent))
% title('PC normalized variance')
% 
% figure
% imagesc(coeff)
% set(gca,'YTick',1:size(invitro_norm,2),'YTickLabels',invitro_vars)
% set(gca,'XTick',1:size(invitro_norm,1))
% set(gca,'TickLabelInterpreter','none')
% xlabel('PCs')
% title('PCA loadings per variable')
%% OFF Cluster the overlap per layer plus soma position
% close all
% 
% %define the number of invitro_clusters
% invitro_clunum = 6;
% %define the PCs to use
% ov_cluster = score(:,1:3);
% % ov_cluster = ov_cluster(:,[2 3]);
% %cluster using linkage and cluster so that I can get the indexes
% l_tree = linkage(ov_cluster,'average','euclidean');
% invitro_clusters = cluster(l_tree,'maxclust',invitro_clunum);
% %code to order the leaves more optimally
% D = pdist(ov_cluster);
% leafOrder = optimalleaforder(l_tree,D);
% %plot the dendrogram
% figure
% dendrogram(l_tree,0,'orientation','left','Reorder',leafOrder)
% 
% %get the limits of the soma center distribution (shifted to not leave the
% %last point at the edge
% soma_lim = [min(soma_cent(:,2))-1,max(soma_cent(:,2))+1];
% %plot maps of a set of invitro_clusters
% %for all the invitro_clusters
% for clu = 1:invitro_clunum
%     %if the cluster only contains 1 cell
%     if sum(invitro_clusters==clu)==1
%         %skip it
%         continue
%     end
%     %average the cells in question
%     clu_cells = squeeze(mean(maps_mat(:,:,invitro_clusters==clu,:)>0,3));
%     
%     %create a figure for the combined plot
%     h = figure;
%     %plot the overlap map
%     subplot(1,8,1:7)
%     map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invitro_clusters==clu))),1,h,10)
%     %get the depths in the cluster
%     soma_clu = soma_cent(invitro_clusters==clu,2);
%     %also plot the distribution of depths in the cluster
%     subplot(1,8,8)
%     errorbar(0,mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'go','MarkerFaceColor','g')
%     hold('on')
% 
%     plot((rand(length(soma_clu),1)*2)-1,soma_clu,'ok')
%     
% %     temp_vec = tar_par(invitro_clusters==clu);
% %     plot(zeros(length(soma_clu(temp_vec))),soma_clu(temp_vec),'or')
% %     plot(zeros(length(soma_clu(~temp_vec))),soma_clu(~temp_vec),'ob')
% 
%     %also plot the functional cells
%     [~,~,ib] = intersect(soma_clu,soma_func(:,2));
%     %create a binary vector with the coordinates of the clu members
%     soma_bin = zeros(length(soma_func),1)==1;
%     soma_bin(ib) = 1;
%     %extract the cells of interest
%     soma_clu_func = soma_func(soma_bin,2);
% %     plot(zeros(length(soma_clu_func)),soma_clu_func,'*g')
%     %plot cells according to the criterion calculated above
%     soma_1 = soma_func(soma_bin&tar_par,2);
%     soma_2 = soma_func(soma_bin&~tar_par,2);
%     plot((rand(length(soma_1),1)*2)-1,soma_1,'*r')
%     plot((rand(length(soma_2),1)*2)-1,soma_2,'*b')
%     set(gca,'YLim',soma_lim,'XLim',[-1.1 1.1])
%     yyaxis right
%     ylabel(strcat('Soma depth,','Functional cells:',num2str(length(soma_clu_func))))
% %     ylabel(strcat('Soma depth,','first cells:',num2str(sum(temp_vec))))
% 
%     
%     set(gca,'YTick',[],'YColor','k')
% 
% end
%% OFF Try GMMs for clustering (using ov_cluster vector from above)

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
%     clu_cell{clu,1} = fitgmdist(ov_cluster,clu_vec(clu),'RegularizationValue',0.0001,...
%         'CovarianceType','diagonal','Replicates',20,'Options',opts);
%     
%     %cluster the data accordingly
%     clu_cell{clu,2} = cluster(clu_cell{clu,1},ov_cluster);
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
% %and the associated number of invitro_clusters
% invitro_clunum = clu_vec(bic_min);
% %get the indexes from the best model
% invitro_clusters = clu_cell{bic_min,2};
% 
% %plot the average maps from these invitro_clusters
% %define the smoothing factor
% sf = 10;
% %get the limits of the soma center distribution (shifted to not leave the
% %last point at the edge
% soma_lim = [min(soma_cent(:,2))-1,max(soma_cent(:,2))+1];
% %plot maps of a set of invitro_clusters
% %for all the invitro_clusters
% for clu = 1:invitro_clunum
%     %if the cluster only contains 1 cell
%     if sum(invitro_clusters==clu)==1
%         %skip it
%         continue
%     end
%     %average the cells in question
%     clu_cells = squeeze(mean(maps_mat(:,:,invitro_clusters==clu,:)>0,3));
%     
%     %create a figure for the combined plot
%     h = figure;
%     %plot the overlap map
%     subplot(1,8,1:7)
%     map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invitro_clusters==clu))),1,h,sf)
%     %get the depths in the cluster
%     soma_clu = soma_cent(invitro_clusters==clu,2);
%     %also plot the distribution of depths in the cluster
%     subplot(1,8,8)
%     plot(zeros(length(soma_clu)),soma_clu,'ok')
%     hold('on')
%     errorbar(0,mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'ro','MarkerFaceColor','r')
%     set(gca,'YLim',soma_lim,'XLim',[-1 1])
%     yyaxis right
%     ylabel('Soma depth')
%     
%     set(gca,'YTick',[],'YColor','k')
% 
% end
%% Calculate horizontal bias
close all
%center the cell maps via interpolation based on the soma information

%define the amount of interpolation (the whole grid is 1035umx1035um, with
%69um between each pair of stimulation points)
int_amount = 70;
%and also the final amplification factor (i.e. size of the map for display)
map_amount = 20;
%allocate memory for the centered maps
invitro_maps = zeros(size(maps_mat).*[map_amount map_amount 1 1]);
%and for the side bias index
side_bias = zeros(cell_num,2);
%load the cell's soma x coordinate
soma_x = soma_cent(:,1);
%for all the cells
for cells = 1:cell_num

    %for exc and inh
    for eiv = 1:2
        %load the map
        temp_map = maps_mat(:,:,cells,eiv);
        %use imresize for the interpolation
        interpol_map = imresize(temp_map,int_amount);
        %then shift the matrix to put the soma center in the center of the
        %image, and erase the portion that shifted around
        cent_curr = circshift(interpol_map,-int16(soma_x(cells)).*70/69,2);
        
        %calculate the side bias index
        left_side = double(cent_curr(:,1:size(cent_curr,2)/2)>0);
        right_side = double(cent_curr(:,1+size(cent_curr,2)/2:end)>0);
        side_bias(cells,eiv) = (sum(right_side(:)) - sum(left_side(:)))/...
            (sum(right_side(:)) + sum(left_side(:)));
        %store the map
        invitro_maps(:,:,cells,eiv) = imresize(cent_curr,map_amount/int_amount);
    end
end

%average the indexes for E and I
side_bias = squeeze(mean(side_bias,2));

%eliminate NaNs to 0
side_bias(isnan(side_bias)) = 0;

%plot the index distributions
figure
% subplot(1,2,1)
histogram(side_bias(:,1))
% subplot(1,2,2)
% histogram(side_bias(:,2))
% %compare excitation vs inhibition side bias on a cell to cell basis
% figure
% subplot(1,2,1)
% %for all the cells
% for cells = 1:cell_num
%     %plot each pair
%     plot([1 2],side_bias(cells,:),'o-k')
%     hold('on')
% end
% plot([1 2],mean(side_bias,1),'or-')
% set(gca,'XLim',[0 3])
% subplot(1,2,2)
% histogram(diff(side_bias,1,2))

%plot an example shifted cell
figure
A = imfuse(interpol_map,cent_curr,'ColorChannels',[1 2 0]);
imagesc(A)
%% Concatenate the side bias with the rest of the indexes

% invitro_norm = normr_2([ov_perlayer,cs_all,soma_norm,side_bias],2);
invitro_raw = [ov_perlayer,cs_all,soma_cent(:,2),side_bias];
%correct the pial distance to normalized pial distance (divide by 1035
%microns, distance to white matter)
invitro_raw(:,6) = invitro_raw(:,6)./1035;
invitro_norm = zscore(invitro_raw);
%remove NaNs
invitro_norm(isnan(invitro_norm)) = 0;
%define the variable names
invitro_vars = {'oL1','oL2/3','oL4','oDL','CS','PD','SB'};

%plot the distribution of each column
figure
%get the number of columns
col_num = size(invitro_norm,2);
%for all the columns
for cols = 1:col_num
    subplot(round(sqrt(col_num)),ceil(sqrt(col_num)),cols)
    histogram(invitro_raw(:,cols),10)
    title(invitro_vars{cols})
end
%% PCA the data

close all
%define the variable names
% invitro_vars = {'OL1','OL2/3','OL4','ODL','CL1','CL2/3','CL4','CDL','PD'};

%run a PCA on the normalized data
[invitro_coeff,invitro_score,invitro_latent] = pca(invitro_norm);

figure
imagesc(invitro_score)
ylabel('Cells')
xlabel('PC space variables')
title('PC-weighted data')

figure
plot(invitro_latent./sum(invitro_latent))
hold('on')
plot(cumsum(invitro_latent./sum(invitro_latent)))
title('PC normalized variance')

figure
imagesc(invitro_coeff)
set(gca,'YTick',1:size(invitro_norm,2),'YTickLabels',invitro_vars)
set(gca,'XTick',1:size(invitro_norm,1))
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
title('PCA loadings per variable')

%plot the first 3 PCs
figure
plot3(invitro_score(:,1),invitro_score(:,2),invitro_score(:,3),'*')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
%% Cluster the cells including the side bias index
close all
%define the number of invitro_clusters
invitro_clunum = 6;

%select the target PCs
ov_cluster = invitro_score(:,1:2);
%make NaNs 0s
ov_cluster(isnan(ov_cluster)) = 0;

% ov_cluster = ov_cluster(1:48,:);

% %add the normalized depth information
% soma_norm = normr_2(soma_cent(:,2));
% side_norm = normr_2(side_bias);
% ov_cluster = cat(2,ov_cluster,soma_norm,side_norm);

% ov_cluster = ov_cluster(:,[2 3]);
%cluster using linkage and cluster so that I can get the indexes
l_tree = linkage(ov_cluster,'average','euclidean');
invitro_clusters = cluster(l_tree,'maxclust',invitro_clunum);
%code to order the leaves more optimally
D = pdist(ov_cluster);
leafOrder = optimalleaforder(l_tree,D);
%plot the dendrogram
figure
dendrogram(l_tree,0,'orientation','top','Reorder',leafOrder)

%get the limits of the soma center distribution (shifted to not leave the
%last point at the edge
soma_lim = [min(soma_cent(:,2))-1,max(soma_cent(:,2))+1];
%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invitro_clunum
    %if the cluster only contains 1 cell
    if sum(invitro_clusters==clu)==1
        %skip it
        continue
    end
    %average the cells in question
%     clu_cells = squeeze(mean(invitro_maps(:,:,invitro_clusters==clu,:)>0,3));
    clu_cells = squeeze(mean(invitro_maps(:,:,invitro_clusters==clu,:),3));
    
    %create a figure for the combined plot
    h = figure;
    %plot the overlap map
%     subplot(16,2,1:2:32)
    subplot(1,2,1)
    map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invitro_clusters==clu))),1,h,1,1)
    %get the depths in the cluster
    soma_clu = soma_cent(invitro_clusters==clu,2);
    %extract the side bias
    side_clu = side_bias(invitro_clusters==clu,:);
    %also plot the distribution of depths in the cluster
%     subplot(16,2,[4 6 8 10])
    subplot(1,2,2)
    plot(side_clu,soma_clu,'ok')
    hold('on')
    
    %and plot
    
    %also plot the functional cells
    [~,ia,ib] = intersect(soma_clu,soma_func(:,2));
    %create a binary vector with the coordinates of the clu members
    soma_bin = zeros(length(soma_func),1)==1;
    soma_bin(ib) = 1;
    %extract the cells of interest
    soma_clu_func = soma_clu(ia);
    %extract the side bias from the tar_par vector (selected above
    %from the functional cells)
    side_clu_func = side_clu(ia);

    %plot cells according to the criterion calculated above
    soma_1 = soma_clu_func(tar_par(ib));
    side_1 = side_clu_func(tar_par(ib));
    soma_2 = soma_clu_func(~tar_par(ib));
    side_2 = side_clu_func(~tar_par(ib));
    plot(side_1,soma_1,'*r')
    plot(side_2,soma_2,'*b')
    set(gca,'YLim',soma_lim,'XLim',[-1.1 1.1],'Ydir','reverse')
    %plot a cross in 0,0
    plot(zeros(2,1),get(gca,'YLim'),'-k')
    plot(get(gca,'XLim'),[sum(soma_lim)/2,sum(soma_lim)/2],'-k')
    %plot the average
    errorbar(mean(side_clu(:,1)),mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'go','MarkerFaceColor','g','MarkerSize',5,'LineWidth',3)
    yyaxis right
    ylabel(strcat('Soma depth,','Functional cells:',num2str(length(soma_clu_func))))
    xlabel('Side bias')
%     ylabel(strcat('Soma depth,','first cells:',num2str(sum(temp_vec))))

    
    set(gca,'YTick',[],'YColor','k')

end
%% Calculate the cluster averages
close all
%allocate memory for the cluster average
invitro_cluave = zeros(invitro_clunum,size(invitro_norm,2));
%and for the number of cluster members
invitro_clumem = zeros(invitro_clunum,1);
%for all the clusters
for clu = 1:invitro_clunum
    %calculate the average of this cluster
    invitro_cluave(clu,:) = mean(invitro_norm(invitro_clusters==clu,:),1);
    %and the number of members
    invitro_clumem(clu) = sum(invitro_clusters==clu);
end

figure
imagesc(invitro_cluave)
%% Save the calculations

%define the save path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_invitro.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(save_path,save_name),'invitro_raw','invitro_norm','invitro_score','invitro_clunum',...
    'invitro_clusters','invitro_coeff','invitro_latent','invitro_maps',...
    'invitro_names','invitro_vars','invitro_cluave','invitro_clumem')
%% OFF Look for correlation against one of the variables in the functional structure
% close all
% 
% %extract the clustering vector from only the functional cells
% [~,ia,ib] = intersect(all_cell{1}.names,all_cell{3}.names);
% func_cluster = ov_cluster(ia,:);
% %define the cluster field names
% clu_names = {'L1Over','L2/3Over','L4Over','DLOver','PiaD','ExcSideBias','InhSideBias'};
% %and the matching ones from the functional one (not sure why there are 5
% %mismatches)
% func_info2 = func_info(ib);
% %get a vector with the fields that have a scalar per cell (no arrays or
% %anything else)
% names = fieldnames(func_info2);
% 
% %allocate memory for the correlations
% % func_corr = zeros(length(names),1);
% func_corr = cell(length(names),2);
% %allocate memory to store the position of the fields
% use_field = zeros(length(names),1);
% 
% %cycle through the fields
% for fieldvar = 1:length(names)
%     %check whether it is a scalar
%     use_field(fieldvar) = isscalar(func_info2(1).(names{fieldvar}));
%     if use_field(fieldvar) == 1
%         %if it is, extract the numbers and get a correlation
%         [func_corr{fieldvar,1},func_corr{fieldvar,2}] = corr(func_cluster,vertcat(func_info2.(names{fieldvar})));
%     end
% end
% 
% %get the number of useful fields times the number of cluster variables
% use_num = sum(use_field)*size(func_cluster,2);
% %transfer only the fields that have something to a matrix
% corr_vec = vertcat(func_corr{use_field==1,1});
% %also the p values
% p_vec = vertcat(func_corr{use_field==1,2});
% %threshold the p-values
% p_plot = corr_vec;
% p_plot(p_vec>0.05) = NaN;
% %plot the results
% figure
% plot(1:use_num,corr_vec,'*b')
% hold('on')
% plot(1:use_num,p_plot,'or')
% 
% %pull out the names of the thresholded variable combinations
% [clu_coord,name_coord] = ind2sub([size(func_cluster,2),sum(use_field)],find(~isnan(p_plot)));
% %build name cells with the involved parameters
% corr_names = names(use_field==1);
% pair_names = [clu_names(clu_coord)',corr_names(name_coord)];
% 
% %use to check for the actual correlation value
% func_info3 = rmfield(func_info2,names(use_field==0));
% figure
% %for the selected fields
% for fields = 1:length(name_coord)
%     %show the correlation between the variables
%     corr_val = corr(vertcat(func_info3.(corr_names{name_coord(fields)})),func_cluster(:,clu_coord(fields)));
%     %also plot the data
%     clu_data = func_cluster(:,clu_coord(fields));
%     func_data = vertcat(func_info3.(corr_names{name_coord(fields)}));
%     subplot(round(sqrt(length(clu_coord))),ceil(sqrt(length(clu_coord))),fields)
%     
%     loglog(clu_data,func_data,'*')
%     
%     xlabel(strcat(pair_names{fields,1},', corr:',num2str(corr_val)),'Interpreter','none')
%     ylabel(pair_names{fields,2},'Interpreter','none')
%     
%     
% end
%% Normalize and PCA the functional data

close all

%extract the clustering vector from only the functional cells
[~,~,ib] = intersect(all_cell{1}.names,all_cell{3}.names);

%and the matching ones from the functional one (not sure why there are 5
%mismatches)
func_info2 = func_info(ib);

%get rid of the non-responsive cells (non-visually responsive at least)
%get the non responsive vector
resp_vec = vertcat(func_info2(:).responder);
func_info2 = func_info2(resp_vec==1);

%get a vector with the fields that have a scalar per cell (no arrays or
%anything else)
names = fieldnames(func_info2);

%allocate memory for the correlations
func_corr = cell(length(names),1);
%allocate memory to store the position of the fields
use_field = zeros(length(names),1);

%cycle through the fields
for fieldvar = 1:length(names)
    %check whether it is a scalar
    use_field(fieldvar) = isscalar(func_info2(1).(names{fieldvar}));
    if use_field(fieldvar) == 1
        %if it is, extract the numbers and get a correlation
        func_corr{fieldvar,1}= vertcat(func_info2.(names{fieldvar}));
    end
end

% %get the number of useful fields times the number of cluster variables
% use_num = sum(use_field)*size(func_cluster,2);
%transfer only the fields that have something to a matrix
invivo_raw = horzcat(func_corr{use_field==1,1});
% invivo_norm = horzcat(func_cluster,invivo_norm);
%do the same with the names
% invivo_vars = cat(1,ov_names',names(use_field==1));
invivo_vars = names(use_field==1);
 
%normalize the matrix per column
% invivo_norm = normr_2(invivo_norm,2);
invivo_norm = zscore(invivo_raw);
%now exclude binary and NaN fields
invivo_raw = invivo_raw(:,sum(isnan(invivo_norm),1)~=size(invivo_norm,1));
invivo_vars = invivo_vars(sum(isnan(invivo_norm),1)~=size(invivo_norm,1));
invivo_norm = invivo_norm(:,sum(isnan(invivo_norm),1)~=size(invivo_norm,1));

%turn the rest of the NaNs to 0
invivo_norm(isnan(invivo_norm)) = 0;
%kill the binary columns
%get a vector with the columns to keep
nonbin_cols = zeros(size(invivo_norm,2),1);
%for all the columns
for cols = 1:size(invivo_norm,2)
    %determine whether the column is non binary
     [~,ib] = unique(invivo_norm(:,cols));
     nonbin_cols(cols) = size(ib,1)>2;
end
%leave only the non-binary columns
invivo_raw = invivo_raw(:,nonbin_cols==1);
invivo_norm = invivo_norm(:,nonbin_cols==1);
invivo_vars = invivo_vars(nonbin_cols==1);
%perform PCA on the variables
[invivo_coeff,invivo_score,invivo_latent] = pca(invivo_norm);

figure
imagesc(invivo_score)
ylabel('Cells')
xlabel('PC space variables')
title('PC-weighted data')

figure
plot(invivo_latent./sum(invivo_latent))
hold('on')
plot(cumsum(invivo_latent./sum(invivo_latent)))
title('PC normalized variance')

figure
imagesc(invivo_coeff)
set(gca,'YTick',1:size(invivo_norm,2),'YTickLabels',invivo_vars)
set(gca,'XTick',1:size(invivo_norm,1))
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
title('PCA loadings per variable')
%% Cluster the functional data
close all
%define how many PCs to use
ov_cluster = invivo_score(:,1:3);

%define the cluster number
invivo_clunum = 3;
%cluster using linkage and cluster so that I can get the indexes
l_tree = linkage(ov_cluster,'average','euclidean');
invivo_clusters = cluster(l_tree,'maxclust',invivo_clunum);
%code to order the leaves more optimally
D = pdist(ov_cluster);
leafOrder = optimalleaforder(l_tree,D);
%plot the dendrogram
figure
dendrogram(l_tree,0,'orientation','left','Reorder',leafOrder)

%get the limits of the soma center distribution (shifted to not leave the
%last point at the edge
soma_lim = [min(soma_cent(:,2))-1,max(soma_cent(:,2))+1];
%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invivo_clunum
    %if the cluster only contains 1 cell
    if sum(invivo_clusters==clu)==1
        %skip it
        continue
    end
    %average the cells in question
%     clu_cells = squeeze(mean(maps_mat(:,:,invitro_clusters==clu,:)>0,3));
    clu_cells = squeeze(mean(invitro_maps(:,:,invivo_clusters==clu,:),3));
    
    %create a figure for the combined plot
    h = figure;
    %plot the overlap map
    subplot(1,8,1:7)
    map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invivo_clusters==clu))),1,h,10)
    %get the depths in the cluster
    soma_clu = soma_cent(invivo_clusters==clu,2);
    %also plot the distribution of depths in the cluster
    subplot(1,8,8)
    errorbar(0,mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'go','MarkerFaceColor','g')
    hold('on')

    plot((rand(length(soma_clu),1)*2)-1,soma_clu,'ok')
    

    set(gca,'YLim',soma_lim,'XLim',[-1.1 1.1],'Ydir','reverse')
    yyaxis right
    ylabel(strcat('Soma depth'))
%     ylabel(strcat('Soma depth,','first cells:',num2str(sum(temp_vec))))

    
    set(gca,'YTick',[],'YColor','k')

end
%% Calculate the cluster averages
close all
%allocate memory for the cluster average
invivo_cluave = zeros(invivo_clunum,size(invivo_norm,2));
%and for the number of cluster members
invivo_clumem = zeros(invivo_clunum,1);
%for all the clusters
for clu = 1:invivo_clunum
    %calculate the average of this cluster
    invivo_cluave(clu,:) = mean(invivo_norm(invivo_clusters==clu,:),1);
    %and the number of members
    invivo_clumem(clu) = sum(invivo_clusters==clu);
end

figure
imagesc(invivo_cluave)
%% Look at interaction between in vivo and in vitro invitro_clusters
close all
%get the names of the in vivo cells
invivo_names = squeeze(all_cell{3}.names);
invivo_names = invivo_names(vertcat(func_info.responder)==1);

%get the common cells between in vivo and in vitro
[~,ia,ib] = intersect(cell_names,invivo_names);
invivo_names = invitro_names(ia);

%get the properties of the common cells
invitro_cells = invitro_raw(ia,:);
% invivo_cells = invivo_raw;

%create a color map based on the number of invivo invitro_clusters
invivo_color = jet(invivo_clunum);

%get the invitro cluster indexes of the iviv cells
invitro_subclusters = invitro_clusters(ia);

%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invitro_clunum
    %if the cluster only contains 1 cell
    if sum(invitro_clusters==clu)==1
        %skip it
        continue
    end
    %average the cells in question
%     clu_cells = squeeze(mean(invitro_maps(:,:,invitro_clusters==clu,:)>0,3));
    clu_cells = squeeze(mean(invitro_maps(:,:,invitro_clusters==clu,:),3));
    
    %create a figure for the combined plot
    h = figure;
    %plot the overlap map
    subplot(1,2,1)
    map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invitro_clusters==clu))),1,h,1,1)
    %get the depths in the cluster
    soma_clu = soma_cent(invitro_clusters==clu,2);
    %extract the side bias
    side_clu = side_bias(invitro_clusters==clu,:);
    %also plot the distribution of depths in the cluster
    subplot(1,2,2)
    %first the responses
    plot(side_clu,soma_clu,'ok')
    hold('on')
    
    %then the average
    errorbar(mean(side_clu(:,1)),mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'go','MarkerFaceColor','g')
    
    %for all the invivo invitro_clusters
    for fclu = 1:invivo_clunum
        %get the cells from this invivo cluster and plot
        soma_funcclu = invitro_cells(invivo_clusters==fclu&invitro_subclusters==clu,6);
        side_funcclu = invitro_cells(invivo_clusters==fclu&invitro_subclusters==clu,7);
        
        plot(side_funcclu,soma_funcclu,'*','MarkerFaceColor',invivo_color(fclu,:))
    end
    set(gca,'YLim',soma_lim,'XLim',[-1.1 1.1],'Ydir','reverse')
    %plot a cross in 0,0
    plot(zeros(2,1),get(gca,'YLim'),'-k')
    plot(get(gca,'XLim'),[sum(soma_lim)/2,sum(soma_lim)/2],'-k')
    yyaxis right
    ylabel(strcat('Soma depth,','Functional cells:',num2str(sum(invitro_clusters==clu))))
    xlabel('Side bias')
    
    set(gca,'YTick',[],'YColor','k')

end
%% Save the invivo calculations

%define the save path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_invivo.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(save_path,save_name),'invivo_raw','invivo_norm','invivo_score','invivo_clunum',...
    'invivo_clusters','invivo_coeff','invivo_latent','invivo_names',...
    'invivo_vars','invivo_cluave','invivo_clumem')