%% Clean up
clearvars
close all
%% Define constants

%define the layers
layer_c = {[1 2],[3 4 5 6],[7 8 9],10:16};
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
%% Generate overlap maps by binarizing the images
close all
%create the matrix with the maps
bin_map = maps_mat>0;

%define the smooth factor
sf = 1;
%generate a blank matrix to fill in the other color channels
blank = ones(size(bin_map,1),size(bin_map,2));
blankw = ones(size(bin_map,1),size(bin_map,2),3);
%calculate the average across cells for both E and I
ave_map = squeeze(mean(bin_map,3));

% %set the excitation map as a red image, smoothing by sf
exc_map = imresize(cat(3,blank,1-normr_2(ave_map(:,:,1)),1-normr_2(ave_map(:,:,1))),sf);
% %and the inhibition map as a blue one
inh_map = imresize(cat(3,1-normr_2(ave_map(:,:,2)),1-normr_2(ave_map(:,:,2)),blank),sf);
% %blend the two images using alpha (and make it double cause default is 8bit
im_ex = double(imfuse(exc_map,inh_map,'method','blend'));

figure

%normalize the image
im_ex = normr_2(im_ex);
%remove the NaNs (if there was an empty channel)
im_ex(isnan(im_ex)) = 0;

%actually plot the map
imagesc(im_ex)
%% Also calculate overlap index based on the bin maps
close all
%allocate memory for the index
ov_index = zeros(cell_num,1);

%for all the cells
for cells = 1:cell_num
    %get the overlap map
    ov_map = squeeze(sum(bin_map(:,:,cells,:),4));
    %quantify the index
    ov_index(cells) = sum(ov_map(:)>1)./sum(ov_map(:)>0);
end

%plot the distribution of indexes
figure
histogram(ov_index)
title('Overlap Index')
%% Calculate overlap per layer

close all hidden

%allocate memory for the indexes per layer
ov_perlayer = zeros(cell_num,layer_num);

%for all the cells
for cells = 1:cell_num
    %get the overlap map
    ov_map = squeeze(sum(bin_map(:,:,cells,:),4));
    %for all the layers
    for layers = 1:layer_num

        %get the layers of interest
        ov_layers = ov_map(layer_c{layers},:);
        %get the index with the corresponding layers
        ov_perlayer(cells, layers) = sum(ov_layers(:)>1)./sum(ov_layers(:)>0);
    end
end

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
%% Cluster the cells by their overlap per layer
close all hidden

%define the number of clusters
clu_num = 4;

%make NaNs 0s
ov_cluster = ov_perlayer;
ov_cluster(isnan(ov_perlayer)) = 0;

ov_cluster = ov_cluster(:,[2 3]);
%cluster using linkage and cluster so that I can get the indexes
l_tree = linkage(ov_cluster,'average','euclidean');
clusters = cluster(l_tree,'maxclust',clu_num);
%code to order the leaves more optimally
D = pdist(ov_cluster);
leafOrder = optimalleaforder(l_tree,D);
figure
dendrogram(l_tree,0,'orientation','left','Reorder',leafOrder)

%plot maps of a set of clusters
%for all the clusters
for clu = 1:clu_num
    %average the cells in question
    clu_cells = squeeze(mean(maps_mat(:,:,clusters==clu,:)>0,3));
    
    %plot the overlap map
    map_plot(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(clusters==clu))),1)
end