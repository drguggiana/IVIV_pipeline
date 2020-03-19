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

%loading path
load_path = 'R:\Share\Simon\Drago_Volker_Simon\';

%file names
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
%% Generate average maps

close all
%generate average images for both types of map
%for both types
for eis = 1:2
    
    mean_map = squeeze(mean(norm_maps(:,:,:,eis),3));
    std_map = squeeze(std(norm_maps(:,:,:,eis),0,3));
    zscore_map = mean_map./std_map;
    
    map_plot(mean_map,'Average')
    
    map_plot(std_map,'STD')
    
%     figure
%     imagesc(zscore_map)
%     colorbar
%     title('Zscore')
    
end
%% Generate overlap maps by binarizing the images
close all
%create the matrix with the maps
bin_map = maps_mat>0;
% bin_map = maps_mat;
% bin_map = norm_maps;

%define the smooth factor
sf = 1;
%generate a blank matrix to fill in the other color channels
blank = ones(size(bin_map,1),size(bin_map,2));
blankw = ones(size(bin_map,1),size(bin_map,2),3);
%calculate the average across cells for both E and I
ave_map = squeeze(mean(bin_map,3));

% %set the excitation map as a red image, smoothing by sf
% exc_map = imresize(cat(3,normr_2(ave_map(:,:,1)),blank,blank),sf);
exc_map = imresize(cat(3,blank,1-normr_2(ave_map(:,:,1)),1-normr_2(ave_map(:,:,1))),sf);
% %and the inhibition map as a blue one
% inh_map = imresize(cat(3,blank,blank,normr_2(ave_map(:,:,2))),sf);
inh_map = imresize(cat(3,1-normr_2(ave_map(:,:,2)),1-normr_2(ave_map(:,:,2)),blank),sf);
% %blend the two images using alpha (and make it double cause default is 8bit
im_ex = double(imfuse(exc_map,inh_map,'method','blend'));

figure

%normalize the image
im_ex = normr_2(im_ex);
%remove the NaNs (if there was an empty channel)
im_ex(isnan(im_ex)) = 0;

% %add a baseline to the image
% im_ex = im_ex + 0.8;

% %invert the image
% im_ex = 1 - im_ex;
% im_ex(:,:,2) = 0;

% %create a mask for the active pixels
% im_mask = reshape(sum(im_ex,3)>0,v_num*h_num*sf*sf,1);
% %reshape the image to invert pixels
% im_ex = reshape(im_ex,v_num*h_num*sf*sf,3);
% %invert only the non-active pixels to get a white background
% im_ex(~im_mask,:) = 1-im_ex(~im_mask,:);
% %reshape the image back
% im_ex = reshape(im_ex,v_num*sf,h_num*sf,3);

imagesc(im_ex)
% cmap = colormap;
% cmap(1,:) = [1 1 1];
% colormap(cmap)
%% Also calculate overlap index based on the bin maps
close all
%allocate memory for the index
ov_index = zeros(cell_num,1);

%for all the cells
for cells = 1:cell_num
%     %get the exc and inh maps
%     exc_map = bin_map(:,:,cells,1);
%     inh_map = bin_map(:,:,cells,2);
%     %calculate the index using the binary maps
%     ov_index(cells) = (sum(exc_map(:)) - sum(inh_map(:)))./...
%         (sum(exc_map(:)) + sum(inh_map(:)));

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
%         %get the current layers for this cell
%         curr_exc = bin_map(layer_c{layers},:,cells,1);
%         curr_inh = bin_map(layer_c{layers},:,cells,2);
%         %calculate the index
%         ov_perlayer(cells,layers) = (sum(curr_exc(:))-sum(curr_inh(:)))./...
%             (sum(curr_exc(:))+sum(curr_inh(:)));
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

% CGobj = clustergram(ov_perlayer,'ImputeFun',@knnimpute,'Cluster','column',...
%     'Colormap','parula','Symmetric',false,'DisplayRange',max(ov_perlayer(:)),...
%     'ColumnLabels',layer_name);

%define the number of clusters
clu_num = 6;

%make NaNs 0s
ov_cluster = ov_perlayer;
ov_cluster(isnan(ov_perlayer)) = 0;
%add the depth information
ov_cluster = cat(2,ov_cluster,soma_cent(:,2));

% ov_cluster = ov_cluster(:,[2 3]);
%cluster using linkage and cluster so that I can get the indexes
l_tree = linkage(ov_cluster,'average','euclidean');
clusters = cluster(l_tree,'maxclust',clu_num);

D = pdist(ov_cluster);
leafOrder = optimalleaforder(l_tree,D);
%plot the dendrogram
figure
dendrogram(l_tree,0,'orientation','left','Reorder',leafOrder)

%plot maps of a set of clusters
%for all the clusters
for clu = 1:clu_num
    %average the cells in question
    clu_cells = squeeze(mean(maps_mat(:,:,clusters==clu,:)>0,3));
    
    %plot the overlap map
    subplot(1,8,1:7)
    map_plot(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(clusters==clu))),1)
    hold('on')
    %get the depths in the cluster
    soma_clu = soma_cent(clusters==clu,2);
    %also plot the distribution of depths in the cluster
    subplot(1,8,8)
    plot(zeros(length(soma_clu)),soma_clu,'*')
end
%% OFF Generate null distribution of overlap for each cell

% %generate null distribution for overlap of signals
% %define the number of iterations
% num_iter = 1000;
% %get the number of points per map
% pts = v_num*h_num;
% %allocate memory for the null results
% ov_null = zeros(v_num*h_num,cell_num);
% %and for the confidence intervals
% ov_conf = zeros(cell_num,1);
% 
% %reshape the norm maps to make space 1D
% oned_maps = reshape(norm_maps,v_num*h_num,cell_num,2);
% %for all the iterations
% for cells = 1:cell_num
%     %show the current cell
%     fprintf(strcat('Current cell:',num2str(cells),'\r\n'))
%     %allocate temp memory for this cell
%     temp_mat = zeros(v_num*h_num,num_iter);
%     %load the current cell
%     curr_cell_e = oned_maps(:,cells,1);
%     curr_cell_i = oned_maps(:,cells,2);
%    
%     %for all the cells
%     for iters = 1:num_iter
%         temp_mat(:,iters) = 1./((curr_cell_e(randperm(pts))-curr_cell_i(randperm(pts))).^2);
%     end
%     temp_mat(isinf(temp_mat(:))) = 0;
%     temp_mat(temp_mat(:)>1000) = 0;
%     %store the average of the iterations
%     ov_null(:,cells) = mean(temp_mat,2);
%     %and calculate the confidence interval
%     ov_conf(cells) = prctile(mean(temp_mat,1),95);
% end
% 
% %return the matrix to its original form
% ov_null = reshape(ov_null,v_num,h_num,cell_num);
% 
% % %plot the average null map
% % map_plot(mean(ov_null,3),'Null overlap map')
%% OFF Calculate overlap between E/I (based on nulls from above)

% %calculate the overlap index for each pair of images
% ov_index = 1./((norm_maps(:,:,:,1)-norm_maps(:,:,:,2)).^2);
% %make the inf values (coming from 0 on both pixels) 0
% ov_index(isinf(ov_index(:))) = 0;
% %also clip values above 1000 since they come from really small differences
% %that explode
% ov_index(ov_index(:)>1000) = 0;
% 
% %show an average overlap index map
% map_plot(mean(ov_index,3),'Overlap Index')
% 
% %clip the overlap based on the average null baseline
% ov_map = mean(ov_index,3);
% ov_map(ov_map<mean(ov_null(:))) = 0;
% map_plot(ov_map,'Clipped baseline overlap')
% 
% %also calculate overlap per cell
% 
% %allocate memory for the overlap
% ov_cells = zeros(cell_num,1);
% %for all the cells
% for cells = 1:cell_num
%     %get the current cell (and 1d it)
%     curr_cell = reshape(ov_index(:,:,cells),v_num*h_num,1);
%     %calculate the overlap excluding values below the threshold for the
%     %cell
%     ov_cells(cells) = mean(curr_cell(curr_cell>ov_conf(cells)));
% end
% 
% %plot the index
% figure
% plot(ov_cells)
% title('Overlap index per cell')
%% OFF Calculate overlap index per layer
% close all
% %allocate memory for the average index per layer per cell
% ov_perlayercell = zeros(cell_num,layer_num);
% 
% %calculate the index
% %for all the cells
% for cells = 1:cell_num
%     %get the current cell
%     curr_cell = ov_index(:,:,cells);
%     
%     %for all the layers
%     for layers = 1:layer_num
%         %get the current layer
%         curr_layer = curr_cell(layer_c{layers},:);
%         %1d it
%         curr_layer = curr_layer(:);
%         %calculate the overlap excluding values below the threshold for the
%         %cell
%         ov_perlayercell(cells,layers) = mean(curr_layer(curr_layer>ov_conf(cells)));
%     end
%     
% end
% % %make NaNs 0 since they come from averaging zero
% % ov_perlayercell(isnan(ov_perlayercell)) = 0;
% 
% %plot a map of cells and layers
% figure
% imagesc(ov_perlayercell')
% set(gca,'YTick',1:4,'YTickLabels',layer_name,'XTick',[])
% xlabel('Cells')
% colorbar
% title('Overlap per layer and cell')
% %also plot the average ov per layer
% figure
% errorbar(1:4,nanmean(ov_perlayercell,1),nanstd(ov_perlayercell,0,1)./sqrt(cell_num),'*')
% title('Overlap per layer')
% set(gca,'XTick',1:4,'XTickLabels',layer_name,'XLim',[0 5])
%% OFF Cluster the cells by their overlap per layer
% close all hidden
% 
% clustergram(ov_perlayercell,'ImputeFun',@knnimpute,'Cluster','column',...
%     'Colormap','parula','Symmetric',false,'DisplayRange',max(ov_perlayercell(:)),...
%     'ColumnLabels',layer_name)
%% Horizontal bias
close all hidden

