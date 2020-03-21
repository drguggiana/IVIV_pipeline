%% Clean up
clearvars
close all
% addpath(genpath('R:\Share\Simon\Drago_Volker_Simon\umap'))
Paths
%% Define the paths

% define the target loading path
% input_path = 'I:\Simon Weiler\INPUT MAPS_final';
input_path = input_maps_path;

% define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\InVivo feature extraction';
% define the path to the morphoMaps
% morpho_path = 'R:\Share\Simon\Drago_Volker_Simon\morphoMaps\morphoMaps.mat';
morpho_path = morphomaps_path;
% define the path to the main data structure
% main_path = 'R:\Share\Simon\Drago_Volker_Simon\InVivivo_InVitro full data structure\190612_1400_dataStruct_ODI_etc_added.mat';
% main_path = 'R:\Share\Simon\Drago_Volker_Simon\_Post_Simon\200219_str_final.mat';
main_path = fullfile(structure_path,'200219_str_final.mat');

% define the path to save the pics
% morpho_pics_path = 'R:\Share\Simon\Drago_Volker_Simon\Morpho_input_pics';

% define the path to the TMD cluster indexes
% tmd_path = 'R:\Share\Simon\Drago_Volker_Simon\TMD_allcells\Analysis_results\cluster_indexes.txt';
tmd_path = tmd_cluster_path;
% define the path to save the aligned maps
aligned_path = 'R:\Share\Simon\Drago_Volker_Simon\InVivo feature extraction\aligned_maps.mat';


%% Load the morpho info
interp_neurons = load(morpho_path);
interp_neurons = interp_neurons.interp_neurons;
%% Load the main data structure

str = load(main_path);
str = str.str;
%% Create the feature vector for the tSNE plot

% define a list of the fields to serialize
% field_list = {'excMap','inhMap','pialD','morph','iv_spon','iv_popcop','iv_OD'};
% field_list = {'excMap','inhMap','pialD','morph','iv_spon','iv_popcop'};
% field_list = {'excMap','inhMap','pialD','iv_spon','iv_popcop'};
% field_list = {'excMap','inhMap','pialD','morph'};
% field_list = {'pialD','morph'};
field_list = {'excMap','inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
iviv_cells = [str(:).iviv]==1;
morpho_cells = ~cellfun(@isempty, {str.morph});
% cell_idx = find(iviv_cells&morpho_cells);
% cell_idx = find(iviv_cells);
% cell_idx = find(morpho_cells);
cell_idx = 1:length(str);
% get the number of cells to include
cell_num = length(cell_idx);

% allocate memory for the individual cells
cell_cell = cell(cell_num,1);


% for all the cells
for cells = 1:cell_num
    % allocate a small cell to hold each field
    field_cell = cell(field_number,1);
    % for all the fields
    for field = 1:field_number
        field_cell{field} = str(cell_idx(cells)).(field_list{field})(:);
        if contains(field_list{field}, 'excMap')
            field_cell{field} = field_cell{field}./min(field_cell{field});
        elseif contains(field_list{field}, 'inhMap')
            field_cell{field} = field_cell{field}./max(field_cell{field});
        end
        
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});
end

% concatenate the results
cell_cell = cat(2,cell_cell{:})';
% remove cells with NaNs
non_nan_cells = sum(isnan(cell_cell),2)==0;
cell_cell = cell_cell(non_nan_cells,:);
% copy the cell to have the original maps later
original_maps = cell_cell;
% redefine cell number based on the rows that didn't contain NaN
cell_num = size(cell_cell,1);
%% Get the morphology density maps

morpho_basal = {str.morphoMap_basal}';
morpho_basal = morpho_basal(non_nan_cells);

morpho_apical = {str.morphoMap_apical}';
morpho_apical = morpho_apical(non_nan_cells);
%% Align the exc and inh maps vertically

% % define the number of bins
% bin_num = 5;
% distance between squares
grid_distance = 69;
% get the square edges
edges = 0:grid_distance:16*grid_distance;
% get the exc and inh maps
% clean_ex = reshape(cat(3,str(non_nan_cells).excMap),256,[]);
% clean_in = reshape(cat(3,str(non_nan_cells).inhMap),256,[]);
clean_ex = cell_cell(:,1:256)';
clean_in = cell_cell(:,257:512)';

% % define the cell number
% cell_num = size(clean_ex,2);
% bin the pial vector
% pia_bins = discretize(cell_cell(:,513),bin_num);
pia_bins = discretize(cell_cell(:,513),edges);

% get the number of actually populated bins
bin_num = max(pia_bins);
% allocate memory for the aligned maps
aligned_maps_ex = zeros(cell_num,16+bin_num,16);
aligned_maps_in = zeros(cell_num,16+bin_num,16);
aligned_maps_basal = zeros(cell_num,16+bin_num,16);
aligned_maps_apical = zeros(cell_num,16+bin_num,16);

% for all the cells
for cells = 1:cell_num
    % get the cells displacement
    cell_soma = pia_bins(cells);   
    % place the map into the aligned matrix based on this displacement
    aligned_maps_ex(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_ex(:,cells),16,16);
    aligned_maps_in(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_in(:,cells),16,16);
    % if the morpho is missing, skip
    if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
        aligned_maps_basal(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
            morpho_basal{cells};        
    end
    % if the morpho is missing, skip
    if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
        aligned_maps_apical(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
            morpho_apical{cells};        
    end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
%% Align the exc and inh maps horizontally

% % define the number of bins
% hbin_num = 5;
% get the exc and inh maps
% clean_ex = reshape(cat(3,str(non_nan_cells).excMap),256,[]);
% clean_in = reshape(cat(3,str(non_nan_cells).inhMap),256,[]);
clean_ex = aligned_maps_ex';
clean_in = aligned_maps_in';
clean_basal = aligned_maps_basal;
clean_apical = aligned_maps_apical;

% get the center of the soma
soma_centers = cat(1,str(:).somaCenter);
soma_centers = soma_centers(non_nan_cells,1)- min(soma_centers(non_nan_cells,1));
% % define the cell number
% cell_num = size(clean_ex,2);
% bin the pial vector
% pia_bins = discretize(cell_cell(:,513),bin_num);
% center_bins = discretize(soma_centers, hbin_num);
center_bins = discretize(soma_centers, edges);

% if all the cells are in the same bin, don't center
if length(unique(center_bins)) == 1
    hbin_num = 0;
else
    % get the number of actually populated bins
    hbin_num = max(center_bins);
    % allocate memory for the aligned maps
    aligned_maps_ex = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_in = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_basal = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_apical = zeros(cell_num,16+bin_num, 16+hbin_num);
    % for all the cells
    for cells = 1:cell_num
        % get the cells displacement
        cell_soma = center_bins(cells);   
        % place the map into the aligned matrix based on this displacement
        aligned_maps_ex(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_ex(:,cells),16+bin_num,16);
        aligned_maps_in(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_in(:,cells),16+bin_num,16);
        % if the morpho is missing, skip
        if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
            aligned_maps_basal(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_basal(cells,:,:);
        end
        % if the morpho is missing, skip
        if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
            aligned_maps_apical(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_apical(cells,:,:);
        end
    end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
%% Plot the aligned maps
close all
figure
plot_count = 1;
for maps = 1:cell_num
%     curr_map = reshape(aligned_maps_ex(maps,:),21,21);
%     curr_map2 = squeeze(aligned_maps_basal(maps,:,:));
%     curr_map = reshape(aligned_maps_ex(maps,:),16+bin_num,16);
%     curr_map2 = squeeze(aligned_maps_basal(maps,:,:)+aligned_maps_apical(maps,:,:));
    curr_map = reshape(aligned_maps_ex(maps,:),16+bin_num,16+hbin_num);
    curr_map2 = squeeze(aligned_maps_basal(maps,:,:)+aligned_maps_apical(maps,:,:));

%     curr_map = reshape(cell_cell(maps,1:256),16,16);
%     curr_map2 = morpho_basal{maps}+morpho_apical{maps};
    if sum(curr_map2(:)) == 0
        continue
    end

    curr_map = cat(3,curr_map,curr_map2,curr_map);
    subplot(5,5,plot_count)
    imagesc(curr_map);
    hold on
    plot([sum(get(gca,'XLim'))/2,sum(get(gca,'XLim'))/2],get(gca,'YLim'),'w')
    axis square
    
    plot_count = plot_count + 1;
    if plot_count > 25
        figure
        plot_count = 1;
    end
end
autoArrangeFigures
% error('stop')
%% Run 2 separate PCAs on ex and in
[coeff_ex,score_ex,latent_ex,~,explained_ex,mu] = pca(aligned_maps_ex);
[coeff_in,score_in,latent_in,~,explained_in,mu] = pca(aligned_maps_in);
%% Save the aligned maps and PCs

% save(aligned_path,'aligned_maps_ex','aligned_maps_in','coeff_ex','coeff_in','score_ex','score_in',...
%     'aligned_maps_basal','aligned_maps_apical')
%% Plot the pca components

% define the number of PCs to plot
pc_num = 3;
figure
% for the first n PCs
for pc = 1:pc_num
    subplot(2,pc_num,pc)
    imagesc(reshape(coeff_ex(:,pc),16+bin_num,16+hbin_num))
    axis square
    subplot(2,pc_num,pc+pc_num)
    imagesc(reshape(coeff_in(:,pc),16+bin_num,16+hbin_num))
    axis square
end

figure
plot(cumsum(explained_ex),'o')
hold on
plot(cumsum(explained_in),'*')
%% Assemble the feature vector

% cell_cell = cat(2,cell_cell(:,513:end),score_ex(:,1:pc_num),score_in(:,1:pc_num));
% cell_cell = cat(2,cell_cell(:,513:end),score_ex(:,1:3),score_in(:,1:3));
 cell_cell = cat(2,score_ex(:,1:3),score_in(:,1:3));
%% Cluster the maps
pialD = [str(:).pialD];
clu_num = 4;
[idx_input, clustering_input, leafOrder] = hca(cell_cell(:,2:end),0,'ward',clu_num,pialD(non_nan_cells),0);%call function for clustering
%% OFF Plot the cluster average maps

% figure
% 
% % for all the clusters
% for clu = 1:clu_num
%     subplot(2,clu_num,clu)
%     imagesc(reshape(mean(original_maps(clu==idx_input,1:256),1),16,16))
%     hold on
%     plot(get(gca,'XLim'),[8.5 8.5],'k')
%     plot([8.5 8.5],get(gca,'YLim'),'k')
%     caxis([0 0.8])
%     axis square
%     subplot(2,clu_num,clu+clu_num)
%     imagesc(reshape(mean(original_maps(clu==idx_input,257:512),1),16,16))
%     hold on
%     plot(get(gca,'XLim'),[8.5 8.5],'k')
%     plot([8.5 8.5],get(gca,'YLim'),'k')
%     caxis([0 0.8])
%     axis square
% end
% 
% figure
% soma_depth = pialD(non_nan_cells);
% 
% % for all the clusters
% for clu = 1:clu_num
%     errorbar(clu,mean(soma_depth(clu==idx_input)),std(soma_depth(clu==idx_input))/sqrt(sum(clu==idx_input)),'o')
%     hold on
% end
% set(gca,'XLim',[0, clu_num+1],'YDir','reverse')
%% Read the TMD cluster indexes

% read the file
idx_table = readtable(tmd_path);
% get the names
tmd_names = table2cell(idx_table(:,1));
% get the matching indexes
tmd_idx = table2array(idx_table(:,2));

% get the cell names from the str
str_names = {str(:).cellName};
str_names = str_names(non_nan_cells);

% build the index vector
% allocate memory for the indexes
clusters = zeros(cell_num,1);
for cells = 1:cell_num
    if sum(contains(tmd_names,str_names{cells})) > 0
        clusters(cells) = tmd_idx(contains(tmd_names,str_names{cells}));
    else
        clusters(cells) = NaN;
    end
end
%% Morpho density corr
% calculate the correlation between morpho density and maps
morphoMaps = {str(non_nan_cells).morphoMap_apical};
excMaps = {str(non_nan_cells).excMap};
inhMaps = {str(non_nan_cells).inhMap};
non_nan_names = {str(non_nan_cells).cellName};

correlation_values = cell(cell_num,3);

% for all the cells
for cells = 1:cell_num
    if isempty(morphoMaps{cells})
        correlation_values{cells,3} = 'nanCell';
        continue
    end
    % save the name of the cell
    correlation_values{cells,3} = non_nan_names{cells};
    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = excMaps{cells}*-1;
            case 2
                map = inhMaps{cells};
        end
        % calculate the correlation and store
        correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
        
    end
    
    

end
%% Save the input features to an output file 

% allocate memory to store the data
input_parameters = struct();

% find the indexes of the non_nan_cells
non_nan_idx = find(non_nan_cells);

% for all the cells
for cells = 1:cell_num
    input_parameters(cells).cellName = str(non_nan_idx(cells)).cellName;
    input_parameters(cells).inputPCs = cell_cell(cells,1:6);
    input_parameters(cells).clusterIdx = idx_input(cells);
end

% save the file
% save(fullfile(out_path,'input_parameters.mat'),'input_parameters')
%% OFF Perform the tsne analysis

% % normalize the columns of the matrix
% cell_cell = normr_2(cell_cell,2);
% % remove the NaNs
% cell_cell(isnan(cell_cell)) = 0;
% 
% tsne_out = tsne(cell_cell);
%% OFF Generate the tSNE plot with a function

% close all
% plot_selector = [1:16];
% var_matrix = cat(2,clusters, idx_input, round(score_ex(:,1)), round(score_in(:,1)));
% var_titles = {'TMD clusters','Input map clusters','First Exc input map PC','First Inh input map PC'};
% plotting_embedding(tsne_out, str, plot_selector,non_nan_idx, non_nan_cells, correlation_values, var_matrix, var_titles)
%% OFF Generate the tsne plot
% % close all
% 
% % % define a color map based on a desired feature
% % parameter = round([str(cell_idx).pialD]);
% % parameter = parameter(non_nan_cells);
% 
% 
% 
% % parameter = round([str(cell_idx).iv_spon].*100);
% 
% % % cross modal parameters
% % iv_param = {str(non_nan_cells).morph};
% iv_param = {str(non_nan_cells).Oripref}';
% % iv_param = correlation_values(:,2);
% parameter = zeros(cell_num,1);
% for cells = 1:cell_num
%     if ~isempty(iv_param{cells}) == 1
%         parameter(cells) = round(iv_param{cells}(1).*100);
%     else
%         parameter(cells) = NaN;
%     end
% end
% 
% % % TMD clusters
% % parameter = clusters;
% 
% % % setup
% % parameter = [zeros(1,45),ones(1,112)];
% % % time
% % parameter = 1:157;
% % % hemisphere
% % parameter = [str(:).hemisphere];
% % % slice orientation
% % parameter = [str(:).sliceOri];
% 
% % parameter = parameter(non_nan_cells);
% 
% 
% % parameter = parameter(non_nan_cells);
% % parameter = round(score_ex(:,1));
% 
% % % input map clusters
% % parameter = idx_input;
% 
% color_number = max(parameter)-min(parameter)+1;
% color_indexer = parameter-min(parameter)+1;
% cmap = parula(color_number);
% 
% figure
% % for all the cells
% for cells = 1:cell_num
%     if ~isnan(parameter(cells))
%         plot(tsne_out(cells,1),tsne_out(cells,2),'o','MarkerFaceColor',cmap(color_indexer(cells),:),'MarkerEdgeColor',cmap(color_indexer(cells),:))
%         hold on
%     else
%         plot(tsne_out(cells,1),tsne_out(cells,2),'*m')
%         hold on
%     end
% end
%% OFF Load the persistence images for the morpho data

% % define the file path
% ph_path = 'R:\Share\Simon\Drago_Volker_Simon\TMD_allcells\Analysis_results\persistence_images.txt';
% % load the contents of the file
% fileID = fopen(ph_path);
% raw_ph = textscan(fileID,'%s');
% raw_ph = raw_ph{1};
% fclose(fileID);
%% OFF Format the ph image data

% % get the number of cells
% morpho_num = size(raw_ph,1);
% 
% % allocate memory for the persistence images
% ph_struct = struct([]);
% 
% % for all the cells
% for cells = 1:morpho_num
%     % read the line from the raw
%     raw_line = strsplit(raw_ph{cells},',');
%     % store the name and the persistence image
%     ph_struct(cells).cellName = raw_line{1};
%     ph_struct(cells).phImage = cell2mat(cellfun(@str2double,raw_line(2:end),'UniformOutput',0));    
% end
%% OFF Plot a persistence image

% % select the target neuron
% target_neuron = 1;
% 
% % get the dimensions of the persistence image
% ph_dimensions = sqrt(size(ph_struct(target_neuron).phImage,2));
% 
% figure
% imagesc(reshape(ph_struct(target_neuron).phImage,ph_dimensions,[]))
%% OFF Load the non-ABI data extraction

% %define the path
% % sw_path = 'R:\Share\Simon\Drago_Volker_Simon\Manuscript2018\Electrophysiology\Ephys structure\ephys.mat';
% sw_path = 'R:\Share\Simon\Drago_Volker_Simon\Manuscript2018\Electrophysiology\Ephys structure\Thesis_structure\ephys_thesis.mat';
% 
% %load the file
% sw_data = load(sw_path);
% sw_data = sw_data.ephys_thesis;
%% OFF Run tSNE on the morpho TMD images

% % normalize the columns of the matrix
% ph_images = normr_2(cat(1,ph_struct.phImage),2);
% 
% % remove the NaNs
% ph_images(isnan(ph_images)) = 0;
% 
% % do PCA on the ph diagrams
% [~,score,latent] = pca(ph_images);
% 
% % run tSNE
% % tsne_morpho_out = tsne(ph_images);
% tsne_morpho_out = tsne(score(:,1:8));
%% OFF Plot the tSNE results

% % load the TMD clustering results
% % read the file
% idx_table = readtable(tmd_path);
% % get the names
% tmd_names = table2cell(idx_table(:,1));
% % get the matching indexes
% tmd_idx = table2array(idx_table(:,2));
% 
% % % get the cell names from the str
% % str_names = {str(:).cellName};
% % str_parameter = {str(:).Oripref};
% % str_parameter = {str(:).pialD};
% % str_parameter = {str(:).Dirpref};
% % str_parameter = {str(:).morph};
% % str_parameter = {str(:).SFTFresp};
% % str_parameter = {str(:).ODresp};
% 
% % % morphology correlation to input and input clusters
% % str_names = correlation_values(:,3);
% % str_parameter = correlation_values(:,1);
% % str_parameter = num2cell(idx_input);
% % str_parameter = num2cell(score_ex(:,1));
% % str_parameter = num2cell(score_in(:,2));
% 
% % ephys parameters
% str_names = sw_data(:,1);
% str_parameter = sw_data(:,4);
% 
% % build the index vector
% % allocate memory for the indexes
% parameter = zeros(morpho_num,1);
% for cells = 1:morpho_num
%     if sum(contains(str_names,ph_struct(cells).cellName)) > 0
%         if isempty(str_parameter{contains(str_names,ph_struct(cells).cellName)})
%             parameter(cells) = NaN;
%         else
%             parameter(cells) = round(str_parameter{contains(str_names,ph_struct(cells).cellName)}(1).*100);
%         end
%     else
%         parameter(cells) = NaN;
%     end
% end
% 
% % % TMD cluster indexes
% % parameter = tmd_idx;
% 
% color_number = nanmax(parameter)-nanmin(parameter)+1;
% color_indexer = parameter-nanmin(parameter)+1;
% cmap = parula(color_number);
% 
% figure
% % for all the cells
% for cells = 1:morpho_num
%     if ~isnan(parameter(cells))
%         plot(tsne_morpho_out(cells,1),tsne_morpho_out(cells,2),'o','MarkerFaceColor',cmap(color_indexer(cells),:),'MarkerEdgeColor',cmap(color_indexer(cells),:))
%         hold on
%     else
%         plot(tsne_morpho_out(cells,1),tsne_morpho_out(cells,2),'*m')
%         hold on
%     end
% end
%% OD decomposition
% for i=1:length(non_nan_idx)
%     if ~isempty(str(non_nan_idx(i)).iv_OD_decom)==1
%         od_decomp(i,:)=str(non_nan_idx(i)).iv_OD_decom;
%     else
%         od_decomp(i,:)=ones(1,24)*NaN;
%     end
% end
% %% SF decomposition
% for i=1:length(non_nan_idx)
%     if ~isempty(str(non_nan_idx(i)).iv_SFTF_decom)==1
%         sftf_decomp(i,:)=str(non_nan_idx(i)).iv_SFTF_decom;
%     else
%         sftf_decomp(i,:)=ones(1,24)*NaN;
%     end
% end
% %% Eye specificty
% for i=1:length(non_nan_idx)
%     if isempty(str(non_nan_idx(i)).ipsi_only)==1 && ~isempty(str(non_nan_idx(i)).contra_only)==1 && ~isempty(str(non_nan_idx(i)).bino)==1
%         eyes(i,:)=NaN;
%     elseif ~isempty(str(non_nan_idx(i)).contra_only)==1 
%        eyes(i,:)=1;
%     elseif ~isempty(str(non_nan_idx(i)).ipsi_only)==1
%         eyes(i,:)=2;
%     elseif ~isempty(str(non_nan_idx(i)).bino)==1
%         eyes(i,:)=3;
%     end
% end
% eyes(find(eyes==0))=NaN;
% %% Peak horizontal offset per layer
% for i=1:length(non_nan_idx)
%     hori_l(i,:)=str(non_nan_idx(i)).hori_peak_pl;
%     hore_l(i,:)=str(non_nan_idx(i)).hore_peak_pl;
%     diffe_l(i,:)=str(non_nan_idx(i)).diff_hori_peak_pl;
%     hori_frac(i,:)=str(non_nan_idx(i)).frac_inh;
%     hore_frac(i,:)=str(non_nan_idx(i)).frac_exh;
% end
% %% Left or right basal morpho density 
% 
% LMO = zeros(length(non_nan_idx),1);
% RMO = zeros(length(non_nan_idx),1);
% for i=1:length(non_nan_idx)
%     
%     if ~isempty(str(non_nan_idx(i)).morphoMap_basal)==1
% %         LMO(i)=sum(sum(str(non_nan_idx(i)).morphoMap_basal(4:6,4:7)));
% %         RMO(i)=sum(sum(str(non_nan_idx(i)).morphoMap_basal(4:6,9:12)));
%         LMO(i) = sum(sum(str(non_nan_idx(i)).morphoMap_basal(:,1:8)));
%         RMO(i) = sum(sum(str(non_nan_idx(i)).morphoMap_basal(:,9:16)));
%     else
%         LMO(i)=NaN;
%         RMO(i)=NaN;
%     end
% end
% rm_lm=(RMO-LMO)./(RMO+LMO);
%% Run UMAP on the data

% [reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 20, 'min_dist', 0.79);
[reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 10, 'min_dist', 0.1);
%% OFF Plot UMAP results with a function

% close all;
% plot_selector = [1:48];
% var_matrix = cat(2,clusters, eyes,hori_l(:,5),hore_l(:,5),diffe_l(:,5),round(od_decomp(:,4).*100),round(od_decomp(:,12).*100)...
%     ,round(sftf_decomp(:,23).*100),round(sftf_decomp(:,24).*100), idx_input, round(score_ex(:,1))...
%     , round(score_in(:,1)),round(score_ex(:,2)), round(score_in(:,2)),round(score_ex(:,3))...
%     , round(score_in(:,3)),round(rm_lm').*100, round(oripref').*100, round(dirpref').*100);
% % var_matrix = cat(2,clusters, eyes,hori_l(:,5),hore_l(:,5),diffe_l(:,5),round(od_decomp(:,4).*100),round(od_decomp(:,12).*100)...
% %     ,round(sftf_decomp(:,4).*100),round(sftf_decomp(:,12).*100), idx_input, round(data_w_input(:,1))...
% %     , round(data_w_input(:,1)),round(data_w_input(:,2)), round(data_w_input(:,2)),round(data_w_input(:,3))...
% %     , round(data_w_input(:,3)));
% var_titles = {'TMD clusters','eye','Hori Inh peak L4','Hori Exc peak L4','Diff ex in peak L4','OD decomp 1','OD decomp 11',...
%     'SFTF decomp 1','SFTF decomp 11','Input map clusters','First Exc input map PC','First Inh input map PC'...
%     ,'second Exc input map PC','second Inh input map PC','third Exc input map PC','third Inh input map PC',...
%     'Difference Apical density R/L','Orientation preference (°)','Direction preference (°)'};
% plotting_embedding(reduced_data, str, plot_selector, non_nan_idx, non_nan_cells, correlation_values, var_matrix, var_titles);
%% OFF Plot the distance matrices for 2 maps

% % image to image correlation ?
% 
% % close all
% 
% figure
% % hist3(reduced_data)
% h = histogram2(reduced_data(:,1),reduced_data(:,2),'DisplayStyle','tile','ShowEmptyBins','on');
% h.NumBins = [10 10];
% axis square
% 
% % [~,sort_vector] = sort(cell_cell(:,7));
% % [~,sort_vector] = sort(score_ex(:,1));
% [~,sort_vector] = sort(reduced_data(:,1));
% 
% umap_distance = -pdist(reduced_data(sort_vector,:))';
% exc1_distance = -pdist(score_ex(sort_vector,2))';
% inh1_distance = -pdist(score_in(sort_vector,1))';
% pial_distance = -pdist(cell_cell(sort_vector,7))';
% 
% 
% figure
% subplot(2,2,1)
% imagesc(squareform(umap_distance))
% title(strcat('UMAP distance'))
% axis square
% 
% subplot(2,2,2)
% imagesc(squareform(exc1_distance))
% axis square
% title(strcat('First Exc PC:',num2str(corr(umap_distance,exc1_distance))))
% 
% subplot(2,2,3)
% imagesc(squareform(inh1_distance))
% axis square
% title(strcat('First Inh PC:',num2str(corr(umap_distance,inh1_distance))))
% 
% subplot(2,2,4)
% imagesc(squareform(pial_distance))
% axis square
% title(strcat('Pial distance:',num2str(corr(umap_distance,pial_distance))))
%% OFF Test Isomap/MDS

% % compute the distance matrix
% D = squareform(pdist(cell_cell));
% options = struct([]);
% options(1).dims = 1:10;
% 
% [Y, R, E] = IsoMap(D, 'k', 10, options);
% reduced_data = Y.coords{2}';
% 
% % Y = cmdscale(D,2);
% % reduced_data = Y;
%% Fraction L4

% get the fractions
frac = cat(1,str.frac_exv);

% select and average layer 4 only
frac = mean(frac(:,6:7),2);

frac = frac(non_nan_cells);

frac = -1*log(frac);

frac(isinf(frac)) = nan;
%% Centroid angle

cent_angle = cat(1,str.ang_wmapil3);
%% Dir pref

dir_pref = {str.Dir};
dir_vector = zeros(length(dir_pref),1);
for el = 1:length(dir_pref)
    if isempty(dir_pref{el})
        dir_vector(el) = nan;
    else
        dir_vector(el) = dir_pref{el}(1);
    end
    
end
dir_vector = dir_vector(non_nan_cells);

dir_cutoff = {str.DSI};
dir_cut = zeros(length(dir_cutoff),1);
for el = 1:length(dir_cutoff)
    if isempty(dir_cutoff{el})
        dir_cut(el) = nan;
    else
        dir_cut(el) = dir_cutoff{el}(1);
    end
    
end
dir_cut = dir_cut(non_nan_cells);

cutoff_dsi = 0.3;

dir_vector(dir_cut<cutoff_dsi) = nan;
%% Ori cutoff

ori_vector = cat(1,str.Oripref);

ori_cutoff = {str.OSI};
ori_cut = zeros(length(ori_cutoff),1);
for el = 1:length(ori_cutoff)
    if isempty(ori_cutoff{el})
        ori_cut(el) = nan;
    else
        ori_cut(el) = ori_cutoff{el}(1);
    end
    
end
ori_cut = ori_cut(non_nan_cells);

cutoff_dsi = 0.3;

ori_vector(ori_cut<cutoff_dsi) = nan;
%% Selection vector
sel_vec = non_nan_cells;
sel_vec = ~isnan(ori_vector);
%% Plot UMAP/Isomap

close all

% plot(Y.coords{2}(1,:), Y.coords{2}(2,:), 'o')
plot_selector = [1:2,5:9, 13:15,17:18, 30:41];

% var_matrix = cat(2, clusters, idx_input, round(100*score_ex(:,1))...
%     , round(100*score_in(:,1)),round(100*score_ex(:,2)),round(100*score_in(:,2)),round(100*score_ex(:,3))...
%     , round(100*score_in(:,3)),round(rm_lm.*100));
% 
% var_titles = {'TMD clusters', 'Input map clusters','First Exc input map PC','First Inh input map PC'...
%     ,'second Exc input map PC','second Inh input map PC','third Exc input map PC','third Inh input map PC',...
%     'Difference Apical density R/L'};


var_matrix = cat(2, clusters, idx_input, round(100*score_ex(:,1))...
    , round(100*score_in(:,1)),round(100*score_ex(:,2)),round(100*score_in(:,2)),round(100*score_ex(:,3))...
    , round(100*score_in(:,3)),round(normr_2(frac)*100),round(normr_2(cent_angle)*100),...
    round(normr_2(dir_vector)*100),round(normr_2(ori_vector)*100));

var_titles = {'TMD clusters', 'Input map clusters','First Exc input map PC','First Inh input map PC'...
    ,'second Exc input map PC','second Inh input map PC','third Exc input map PC','third Inh input map PC',...
    'Layer 4 ex fraction','ang_wampil3','dir pref','ori_pref'};


plotting_embedding2(reduced_data, str, plot_selector, non_nan_cells, correlation_values,2, var_matrix, var_titles)

% define the number of PCs to plot
pc_num = 3;
figure
% for the first n PCs
for pc = 1:pc_num
    subplot(2,pc_num,pc)
    imagesc(reshape(coeff_ex(:,pc),16+bin_num,16+hbin_num))
    axis square
    subplot(2,pc_num,pc+pc_num)
    imagesc(reshape(coeff_in(:,pc),16+bin_num,16+hbin_num))
    axis square
end
autoArrangeFigures
%% Arrow plots
error('Stop here')
close all
% plot arrows in each point, pointing to the neighbor in the parameter
% sorted
% define the arrow length
arrow_length = 0.1;
% define a color map for the angles
angle_map = jet(360);

% get the coordinate of the neighbor for each point
% define the target parameter
% tar_param = [str.Oripref]';
% tar_param = pialD(non_nan_cells);

% define the parameters as a cell array
tar_cell = {score_ex(:,1); score_in(:,1); score_in(:,2); pialD(non_nan_cells)'};

% allocate memory to store the results
iso_analysis = struct([]);

% for all the parameters
for tar_count = 1:size(tar_cell,1)
    
    tar_param = tar_cell{tar_count};
    % get the indexes of the sorted values
    [~, order] = sort(tar_param);

    % get the 2d distance between points
    point_distance = squareform(pdist(reduced_data));
    % define the number of neighbors to consider for the averaging
    n_neighbors = 10;
    % define the vector of sorted indexes to compare
    sort_vector = 1:size(reduced_data,1)-1;

    % allocate memory for the angles
    angle_vector = nan(size(reduced_data,1),1);

    % plot the embedding
    figure
    plot(reduced_data(:, 1), reduced_data(:,2),'ok')
    hold on
    % for all but the last point
    for points = 1:size(reduced_data,1)-1
        % get the base coordinates
        base = reduced_data(order==points,:);
        % calculate the angle to the neighbor
        centered = reduced_data(order==(points+1),:) - base;
        angle = atan2(centered(2),centered(1));
        % calculate the sum of the angle to the next n neighbors
        % initialize the counter
        centered = [0,0];
%         % get the neighbors based on the closest index in the sorted vector
%         [~,neighbors] = sort(abs(points-sort_vector));
%         neighbors = neighbors(2:n_neighbors+1);
%         % for all the n neighbors
%         for neighbor = neighbors
%             centered = centered + reduced_data(order==neighbor,:) - base;
%         end
%         % wrap the angle
%         angle = atan2(centered(2),centered(1));
        % calculate the tip based on the defined length
        tip = base + [cos(angle)*arrow_length, sin(angle)*arrow_length];
        start = base + [-cos(angle)*arrow_length, -sin(angle)*arrow_length];
        % calculate the arrow color
        raw_angle = round(unwrap(rad2deg(angle)));
        if raw_angle <= 0
            raw_angle = raw_angle + 360;
        end
        % save the angle
        angle_vector(order==points) = raw_angle;
        % draw the arrow
    %     annotation('arrow',[base(1), tip(1)],[base(2), tip(2)])
        arrow(start, tip,[],[],[],0.1,'facecolor',angle_map(raw_angle,:));
    end

%     figure
%     plot(angle_vector)
    
    % store the angles
    iso_analysis(tar_count).angles = angle_vector;
    %% Calculate the distance between neighbors

    % allocate memory for the distance
    distance_vector = nan(size(reduced_data,1),1);
    % for all but the last point
    for points = 1:size(reduced_data,1)-1
         % get the base coordinates
        base = reduced_data(order==points,:);
        % calculate the distance to the neighbor
        distance = pdist([reduced_data(order==(points+1),:); base]);
        % save the distance
        distance_vector(order==points) = distance;

    end
    % store the distances
    iso_analysis(tar_count).distance = distance_vector;
end
%% Plot the arrow plot results
close all

% define the number of color levels to use
color_levels = 150;
% calculate a delta between properties
delta_prop = iso_analysis(1).distance - iso_analysis(3).distance;
delta_norm = round(normr_2(delta_prop).*color_levels + 1);
% create a colormap for the values
cmap_range = [min(delta_norm), max(delta_norm)];
cmap_number = diff(cmap_range);
cmap = jet(cmap_number+1);

figure
for points = 1:size(reduced_data,1)
    if isnan(delta_norm(points))
        continue
    end
    scatter(reduced_data(points, 1), reduced_data(points,2),[],cmap(delta_norm(points), :), 'filled')
    hold on
end
colormap(jet)
colorbar

figure
histogram(delta_prop)

% get the combinations of parameters
combo_matrix = nchoosek(1:size(tar_cell,1),2);
combo_number = size(combo_matrix,1);
figure
for combos = 1:combo_number
    subplot(round(sqrt(combo_number)),ceil(sqrt(combo_number)),combos)
    plot(iso_analysis(combo_matrix(combos,1)).angles, iso_analysis(combo_matrix(combos,2)).angles, 'o')
    hold on
    plot(normr_2(tar_cell{combo_matrix(combos,1)})*360,normr_2(tar_cell{combo_matrix(combos,2)})*360,'x')
end

% figure
% for param = 1:size(tar_cell,1)
%     subplot(round(sqrt(size(tar_cell,1))),ceil(sqrt(size(tar_cell,1))),param)
%     plot(iso_analysis(param).angles, iso_analysis(param).distance,'o')
%     hold on
% end

figure
for combos = 1:combo_number
    [sorted_param,sorted_idx] = sort(tar_cell{combo_matrix(combos,1)});
    subplot(round(sqrt(combo_number)),ceil(sqrt(combo_number)),combos)
%     plot(iso_analysis(combo_matrix(combos,1)).angles(sorted_idx), 'o')
    polarplot(deg2rad(iso_analysis(combo_matrix(combos,1)).angles(sorted_idx)),...
        iso_analysis(combo_matrix(combos,1)).distance(sorted_idx),'o')
    
    hold on
    polarplot(deg2rad(iso_analysis(combo_matrix(combos,2)).angles(sorted_idx)),...
        iso_analysis(combo_matrix(combos,2)).distance(sorted_idx),'o')
%     plot(iso_analysis(combo_matrix(combos,2)).angles(sorted_idx), 'o')
end
%% Calculate rolling correlation
close all
% define the width of the window
window_width = 20;
% define the number of reps for the randomization
rep_number = 100;
% turn the iv parameters into vectors for plotting
ori_pref = paracell_to_vector({str(non_nan_cells).Oripref}', 1);
dir_pref = paracell_to_vector({str(non_nan_cells).Dirpref}', 1);

% define the parameters as a cell array
tar_cell = {ori_pref; score_ex(:,1); score_in(:,1); score_in(:,2); pialD(non_nan_cells)';dir_pref};
tar_names = {'Oripref','exc_first', 'inh_first', 'inh_second', 'pialD','dir_pref'};

rolling_correlation(tar_cell,tar_names,reduced_data,window_width,rep_number)
tilefigs
%% Save the UMAP results

% save(fullfile(out_path,'umap_out.mat'),'reduced_data','umap')
%% OFF rolling correlation
% close all
% % define the width of the window
% window_width = 20;
% % define the number of reps for the randomization
% rep_number = 100;
% % turn the iv parameters into vectors for plotting
% ori_pref = paracell_to_vector({str(non_nan_cells).Oripref}', 1);
% dir_pref = paracell_to_vector({str(non_nan_cells).Dirpref}', 1);
% 
% % define the parameters as a cell array
% tar_cell = {ori_pref; score_ex(:,1); score_in(:,1); score_in(:,2); pialD(non_nan_cells)';dir_pref};
% tar_names = {'Oripref','exc_first', 'inh_first', 'inh_second', 'pialD','dir_pref'};
% % get the number of parameters
% param_num = size(tar_cell,1);
% % get the combinations of parameters
% combo_matrix = nchoosek(1:param_num,2);
% combo_number = size(combo_matrix,1);
% % allocate memory to hold the combo titles
% combo_title = cell(combo_number,1);
% % create the figures
% fig_line = figure();
% fig_map = figure();
% 
% % for all the parameter combinations
% for combo = 1:combo_number
%     % get the 2 parameters
%     param_1 = tar_cell{combo_matrix(combo,1)};
%     param_2 = tar_cell{combo_matrix(combo,2)};
%     % sort the first parameter and apply the sorting to the second
%     [param_1,idx_sort] = sort(param_1);
%     param_2 = param_2(idx_sort);
%     % get the number of elements
%     num_elements = size(param_1,1);
%     
%     % calculate the moving correlation
%     rolling_corr = (movcorr(param_1,param_2,window_width, 'omitnan'));
%     
%     % calculate the random version
%     % allocate memory for the randomization
%     random_roll = zeros(num_elements,rep_number);
%     % for all the reps
%     for reps = 1:rep_number
%         random_roll(:,reps) = (movcorr(param_1,param_2(randperm(num_elements)),window_width, 'omitnan'));
%     end
% %     % allocate memory for the correlation
% %     rolling_corr = zeros(num_elements,2);
% %     % for all the elements in the parameter
% %     for el = 1:num_elements
% %         % if the elements are at the beginning
% %         if el < window_width
% %             rolling_corr(el,1) = corr(param_1(1:window_width),param_2(1:window_width));
% %         elseif el > window_width && el < num_elements - window_width
% %             rolling_corr(el,1) = corr(param_1(el-window_width:el+window_width),param_2(el-window_width:el+window_width));
% %         else
% %             rolling_corr(el,1) = corr(param_1(end-window_width:end),param_2(end-window_width:end));
% %         end
% %     end
%     % store the combo title
%     combo_title{combo} = strcat(tar_names{combo_matrix(combo,1)},'+',tar_names{combo_matrix(combo,2)});
%     % plot the results
%     figure(fig_line)
%     subplot(round(sqrt(combo_number)),ceil(sqrt(combo_number)),combo)
% %     x_vector = 1:length(param_1);
%     x_vector = param_1;
%     plot(x_vector,rolling_corr,'o')
%     hold on
%     % calculate random_roll sorting
%     [random_mean,random_idx] = sort(mean(random_roll,2));
%     random_std = std(random_roll,0,2);
%     random_std = random_std(random_idx);
%     shadedErrorBar(x_vector,random_mean,random_std)
%     title(combo_title{combo}, 'Interpreter','None')
%     plot(get(gca,'XLim'),[0 0],'k--')
%     figure(fig_map)
%     % get the significantly correlated
%     std_roll = std(random_roll,0,2);
%     sig_correlated = abs(rolling_corr) > abs(std_roll);
%     subplot(round(sqrt(combo_number)),ceil(sqrt(combo_number)),combo)
%     plot(reduced_data(idx_sort,1),reduced_data(idx_sort,2),'ko')
%     hold on
%     plot(reduced_data(idx_sort(sig_correlated),1),reduced_data(idx_sort(sig_correlated),2),'ro','MarkerFaceColor','r')
%     title(combo_title{combo}, 'Interpreter','None')
% 
% end
%% OFF try a quiver plot
% % allocate memory for the position and direction components of the arrow
% quiver_vectors = zeros(size(reduced_data,1)-1, 4);
% % for all the points
% for points = 1:size(reduced_data,1)-1
%     % load the location of the arrow
%     base = reduced_data(order==points,:);
%     % save it for use in the quiver
%     quiver_vectors(order==points, [1 2]) = base;
%     % calculate the angle to the neighbor
%     centered = reduced_data(order==(points+1),:) - base;    
%     angle = atan2(centered(2),centered(1));
%     % calculate the direction components
%     quiver_vectors(order==points, [3 4]) = [cos(angle)*arrow_length, sin(angle)*arrow_length];
% end
% 
% % plot the quiver
% quiver(quiver_vectors(:,1),quiver_vectors(:,2),quiver_vectors(:,3),quiver_vectors(:,4))
%% OFF Plot UMAP results

% % define a color map based on a desired feature
% parameter = round([str(cell_idx).pialD]);
% parameter = parameter(non_nan_cells);

% cross modal parameters
% % iv_param = {str(non_nan_cells).morph};
% iv_param = {str(non_nan_cells).Oripref}';
% % iv_param = correlation_values(:,1);
% parameter = zeros(cell_num,1);
% for cells = 1:cell_num
%     if ~isempty(iv_param{cells}) == 1
%         parameter(cells) = round(iv_param{cells}(1).*100);
%     else
%         parameter(cells) = NaN;
%     end
% end

% % TMD clusters
% parameter = clusters;

% % setup
% parameter = [zeros(1,45),ones(1,112)];
% % time
% parameter = 1:157;
% % hemisphere
% parameter = [str(:).hemisphere];
% % slice orientation
% parameter = [str(:).sliceOri];

% parameter = parameter(non_nan_cells);


% % parameter = parameter(non_nan_cells);
% parameter = round(score_ex(:,1));
% 
% % % input map clusters
% % parameter = idx_input;
% 
% color_number = max(parameter)-min(parameter)+1;
% color_indexer = parameter-min(parameter)+1;
% cmap = parula(color_number);
% 
% figure
% % for all the cells
% for cells = 1:cell_num
%     if ~isnan(parameter(cells))
%         plot(reduced_data(cells,1),reduced_data(cells,2),'o','MarkerFaceColor',cmap(color_indexer(cells),:),'MarkerEdgeColor',cmap(color_indexer(cells),:))
%         hold on
%     else
%         plot(reduced_data(cells,1),reduced_data(cells,2),'*m')
%         hold on
%     end
% end