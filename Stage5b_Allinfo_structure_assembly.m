%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = stage4_full_structure_path;
str = load(find_newer_file(main_path));
str = str.str_allcells;

% get the number of cells
cell_num = length(str);
%% Correct the maps
close all

% plot the maps
% initialize a plot counter
plot_count = 1;
figure
% for all cells
for cells = 1:cell_num
    if isempty(str(cells).excMap)
        continue
    end
    subplot(5,5,plot_count)
    imagesc(str(cells).excMap)
    axis square
    set(gca,'XTick',[],'YTick',[])
    % if the max plots are reached, create a figure and update the counter
    if plot_count == 25
        figure
        plot_count = 0;
    end
    
    % update the plot counter
    plot_count = plot_count + 1;
end

% define the grid spacing in microns
grid_spacing = 69;

% get the map size in microns
map_size = round(16.*grid_spacing);
% get the map limits
map_lim = map_size/2-grid_spacing/2;

% create the grid
[Y,X] = ndgrid(-map_lim:grid_spacing:map_lim,...
    -map_lim:grid_spacing:map_lim);

% get the slice ori cloud centers for setup 2
centroid_vector = zeros(cell_num,1);
% get the map centroids
for cells = 1:cell_num
    
    % get either of the maps. If not available, skip
    if ~isempty(str(cells).excMap)
        % get the map
        map = str(cells).excMap(3:5,:);
    elseif ~isempty(str(cells).inhMap)
        map = str(cells).inhMap(3:5,:);
    else
        continue
    end

    % get the components
    cc = bwconncomp(map);
    rp = regionprops(cc,map,{'Area','WeightedCentroid'});
    % if there's more than one, leave the largest
    if cc.NumObjects > 1
        areas = cat(1,rp.Area);
        [~,idx] = max(areas);
        rp = rp(idx);
    end
    % store the centroid
    centroid_vector(cells) = rp.WeightedCentroid(:,1);
end
% get the sliceori
slice_ori = {str.sliceOri};
slice_ori2 = zeros(length(slice_ori),1);
for el = 1:numel(slice_ori)
    if isempty(slice_ori{el})
        slice_ori2(el) = NaN;
    else
        slice_ori2(el) = slice_ori{el};
    end
end
slice_ori = slice_ori2;
% get the setup
setup = (1:cell_num)';
setup = setup>47;
% get the mean centers
soma_slice = [nanmean(centroid_vector(slice_ori==0&setup==1)),...
    nanmean(centroid_vector(slice_ori==1&setup==1))];

% for all the cells
for cells = 1:cell_num
%     % get the soma centers in x
%     soma = str(cells).somaCenter(1);
%     % flip the sign if it's a sliceOri 0 map
%     if str(cells).sliceOri == 0
%         soma = -soma;
%     end
%     % generate the centered grid
%     center_X = X - soma; 
    if isempty(str(cells).excMap) && isempty(str(cells).inhMap)
        continue
    end
    % select the appropriate offset
    if slice_ori(cells) == 0
        offset = (soma_slice(1)-8.5)*grid_spacing;
        % add the corrected soma position
        str(cells).subpixel_soma = [-str(cells).somaCenter(1),...
            str(cells).somaCenter(2)];
    elseif slice_ori(cells) == 1
        offset = (soma_slice(2)-8.5)*grid_spacing;
        % add the corrected soma position
        str(cells).subpixel_soma = [str(cells).somaCenter(1),...
            str(cells).somaCenter(2)];
    end
    % subtract the offset
    center_X = X + offset;

    % for both maps
    for maps = 1:2
        switch maps
            case 1
                map_type = 'excMap';
            case 2
                map_type = 'inhMap';
        end
        % if it's a setup 1 cell, don't interpolate it
        if setup == 0
            str(cells).(strcat('subpixel_',map_type)) = str(cells).(map_type);
        else
            % get the corresponding map
            map = str(cells).(strcat(map_type));
            % interpolate the map
            interpolant = griddedInterpolant(Y,X,map);
            % replace the map
            new_map = interpolant(Y,center_X);
            % also normalize it
            if maps == 1
                factor = min(new_map(:));
            else
                factor = max(new_map(:));
            end
            % save the normalized map
            str(cells).(strcat('subpixel_',map_type)) = new_map./factor;
            % save the raw map also for examples
            str(cells).(strcat('subpixel_raw_',map_type)) = new_map;
        end
    end
end

% plot centroid histograms
figure
histogram(centroid_vector(slice_ori==0&setup==0))
hold on
histogram(centroid_vector(slice_ori==1&setup==0))

% plot the maps
% initialize a plot counter
plot_count = 1;
figure
% for all cells
for cells = 1:cell_num
    subplot(5,5,plot_count)
    if isempty(str(cells).excMap)
        continue
    end
    imagesc(str(cells).subpixel_excMap)
    axis square
    set(gca,'XTick',[],'YTick',[])
    % if the max plots are reached, create a figure and update the counter
    if plot_count == 25
        figure
        plot_count = 0;
    end
    
    % update the plot counter
    plot_count = plot_count + 1;
end
autoArrangeFigures
%% Add the fraction fields
[frac_exh,~,frac_inh,~,frac_exv,~,frac_inv] = iviv_profiles(1:cell_num,str);
frac_exv_m=[zeros(cell_num,2),frac_exv];
% abs_exv_m=[zeros(cell_num,2),abs_exv];
% layer_assign=[zeros(cell_num,2),(1:cell_num)'];
frac_v=[frac_exv_m,frac_inv];
% frac_v = [frac_exv,frac_inv];
frac_h=[frac_exh,frac_inh];

% add the fields to the structure
% for all the cells
for cells = 1:cell_num
    str(cells).frac_vert = frac_v(cells,:);
    str(cells).frac_horz = frac_h(cells,:);
end
%% Add the span fields

% calculate the span
[ex_spanhL23,ex_spanhL4,ex_spanhL5,in_spanhL23,in_spanhL4,in_spanhL5] = ...
    span_perLayer(str);

% add the span fields
% for all the cells
for cells = 1:cell_num
    str(cells).span = [ex_spanhL23(cells),ex_spanhL4(cells),ex_spanhL5(cells),...
        in_spanhL23(cells),in_spanhL4(cells),in_spanhL5(cells)];
end
%% Re-PCA and re-cluster
close all

% select only the input maps
inputmap_idx = ~cellfun(@isempty,{str.excMap})&~cellfun(@isempty,{str.inhMap})...
    &~cellfun(@(x) all(all(isnan(x))),{str.excMap})&~cellfun(@(x) all(all(isnan(x))),{str.inhMap});
inputmap_num = sum(inputmap_idx);
inputmap_ref = find(inputmap_idx);

% get the maps
exc_maps = reshape(cat(3,str(inputmap_idx).subpixel_excMap),256,inputmap_num);
inh_maps = reshape(cat(3,str(inputmap_idx).subpixel_inhMap),256,inputmap_num);

% align them vertically
% distance between squares
grid_distance = 69;
% get the square edges
edges = 0:grid_distance:16*grid_distance;
% get the pial bins
pia_bins = discretize(cat(1,str(inputmap_idx).pialD),edges);

% get the number of actually populated bins
bin_num = max(pia_bins);
% allocate memory for the aligned maps
aligned_maps_ex = zeros(inputmap_num,16+bin_num,16);
aligned_maps_in = zeros(inputmap_num,16+bin_num,16);

% for all the cells
for cells = 1:inputmap_num
    % get the cells displacement
    cell_soma = pia_bins(cells);   
    % place the map into the aligned matrix based on this displacement
    aligned_maps_ex(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(exc_maps(:,cells),16,16);
    aligned_maps_in(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(inh_maps(:,cells),16,16);
end

% run the pca
[coeff_ex,score_ex] = pca(reshape(aligned_maps_ex,inputmap_num,[]));
[coeff_in,score_in] = pca(reshape(aligned_maps_in,inputmap_num,[]));

% % save the PCs and coeffs in a file
% save(pca_path,'coeff_ex','coeff_in','score_ex','score_in')

% define the clustering input
cluster_input = [score_ex(:,1:3),score_in(:,1:3)];
clu_num = 4;

%including the 6 PCs = pial depth

[idx, clustering, leafOrder] = hca(cluster_input, 0, 'ward',...
    clu_num, cat(1,str(inputmap_idx).pialD), 0, 0.6);

% for all the cells
for cells = 1:cell_num
    if inputmap_idx(cells) == 1
        str(cells).Cluster_id = idx(inputmap_ref==cells);
        str(cells).PCs = cluster_input(inputmap_ref==cells,:);
    else
        str(cells).Cluster_id = NaN;
        str(cells).PCs = NaN;
    end
end
% plot the PCs
figure
% for the pcs
for pc = 1:3
    subplot(2,3,pc)
    imagesc(reshape(coeff_ex(:,pc),16+bin_num,16))
    subplot(2,3,pc+3)
    imagesc(reshape(coeff_in(:,pc),16+bin_num,16))
end

% plot the clusters
figure
% for all the clusters
for clu = 1:clu_num
    subplot(2,clu_num,clu)
    imagesc(squeeze(mean(cat(3,str([str.Cluster_id]==clu).subpixel_excMap),3)))
    hold on
    plot([8.5 8.5],get(gca,'YLim'),'k')
    axis square
    set(gca,'XTick',[],'YTick',[])
    subplot(2,clu_num,clu+clu_num)
    imagesc(squeeze(mean(cat(3,str([str.Cluster_id]==clu).subpixel_inhMap),3)))
    hold on
    plot([8.5 8.5],get(gca,'YLim'),'k')
    axis square
    set(gca,'XTick',[],'YTick',[])
end

figure
pialD = cat(1,str(inputmap_idx).pialD);
% for all the clusters
for clu = 1:clu_num
    errorbar(clu,mean(pialD(idx==clu)),...
        std(pialD(idx==clu))./sqrt(sum(clu==idx)),'o')
    hold on
end

set(gca,'XLim',[0 clu_num+1],'YDir','reverse')
%% Update the centroid angles

% get the pial depth
pialD = cat(1,str(inputmap_idx).pialD);
% get the soma x position
soma_coord = cat(1,str(inputmap_idx).subpixel_soma);
% somax = somax(:,1);
% define the index
% idx_centroid = 1:cell_num;

% for ex and in
for exin = 1:2
    % select the corresponding field
    switch exin
        case 1
            map_type = 'subpixel_excMap';
            field1 = 'ex';
        case 2
            map_type = 'subpixel_inhMap';
            field1 = 'in';
    end
    % get the corresponding maps
    maps = cat(3,str(inputmap_idx).(map_type));
    
    % for both layer groups
    for layer = 1:3
        % select the rows for layer 2/3 and layer 4 respectively
        switch layer
            case 1
                layers = 2:6;
                field2 = 'L23';
                row_shift = 2;
            case 2
                layers = 5:8;
                field2 = 'L4';
                row_shift = 5;
            case 3
                layers = 7:11;
                field2 = 'L5';
                row_shift = 7;
        end

        % get the layer
        map_layers = maps(layers,:,:);
        out_ang = centroid_map(map_layers,soma_coord(:,1),soma_coord(:,2),pialD,row_shift,1);
        out_ang_nonwe = centroid_map(map_layers,soma_coord(:,1),soma_coord(:,2),pialD,row_shift,0);

        % for all the cells
        for cells = 1:cell_num
            if inputmap_idx(cells) == 1
                str(cells).(strcat('ang_',field1,field2)) = out_ang(inputmap_ref==cells,:);
                str(cells).(strcat('ang_',field1,field2,'_nonweighted')) = out_ang_nonwe(inputmap_ref==cells,:);
            end
        end
    end
end
%% Add the noise correlations

% get the list of folders
folder_list_raw = dir(od_svd_path);
%get the number of folders
exp_num = length(folder_list_raw);
% allocate memory for the full paths
folder_list = cell(exp_num,1);
% for all the folders
for folders = 1:exp_num
    folder_list{folders} = fullfile(folder_list_raw(folders).folder,...
        folder_list_raw(folders).name);
end
%allocate memory to store the ROI and neuropil data
svd_OD = cell(exp_num,1);
% also for the noise correlations
noise_OD = cell(exp_num,1);

% allocate memory to also save the number of rois per experiment along with
% the experiment name
cell_number = cell(exp_num,2);

%for all the folders
for experiment = 1:exp_num
    % load the cell containing the hosvd decompositions
    svd_OD{experiment} = load(folder_list{experiment});
    noise_OD{experiment} = svd_OD{experiment}.noise_matrix;
    % collapse the second cell dimension
    noise_OD{experiment} = squeeze(mean(noise_OD{experiment},2));
    svd_OD{experiment} = svd_OD{experiment}.cell_cell;
    
    [~,cell_number{experiment,1}] = fileparts(folder_list{experiment});
    cell_number{experiment,1} = cell_number{experiment,1}(4:8);
    cell_number{experiment,2} = size(svd_OD{experiment},1);
end

% concatenate the data
svd_OD = cat(1,svd_OD{:});
noise_OD = cat(1,noise_OD{:});
% % save the noise correlations to file
% save(noise_path,'noise_OD')
%% Match the iviv cells with the all invivo cells

% load the excel spreadsheet with the matching
[~,~,matching_raw] = xlsread(matching_path);

% get the animal names and their corresponding OD names
matching_names = matching_raw(2:end,[4,9,11]);
% % turn the experiment strings into numbers
matching_names(:,2) = cellfun(@eval, matching_names(:,2),'UniformOutput',0);
% for all the rows (cause NaNs in the cellID field -_-
for rows = 1:size(matching_names,1)
    if ~ischar(matching_names{rows,3})
        matching_names{rows,3} = num2str(matching_names{rows,3});
    end
    matching_names{rows,3} = eval(matching_names{rows,3});
end

% allocate memory to store the pial depth
iviv_vector = cell(exp_num,1);

% for all the experiments
for experiment = 1:exp_num
    % get the current experiment name
    [~,current] = fileparts(folder_list{experiment});
    % find the matching animal name from the iviv names
    iviv_idx = contains(matching_names(:,1),current(4:8));
    % assemble a vector with the cell number
    iviv_temp = zeros(cell_number{contains(cell_number(:,1),current(4:8)),2},2);
    % turn the iviv indexes into ones
    iviv_temp(matching_names{iviv_idx,2},1) = 1;
    % and the cellID fields into their corresponding number (only if not
    % NaN)
    if ~isnan(matching_names{iviv_idx,3})
        iviv_temp(matching_names{iviv_idx,2},2) = matching_names{iviv_idx,3};
    end
    % store the resulting vector in the main cell
    iviv_vector{experiment} = iviv_temp;
end

% concatenate the cell
iviv_vector = cat(1,iviv_vector{:});
% % get the corrected cells only
% iviv_vector = iviv_vector(correct_cells,:);
% iviv_vector = iviv_vector(correct_cells_SFTF,:);
%% Save the noise correlations

% get the noise correlations only for the matched cells
% noise_matched = noise_OD(correct_cells,:,:);
% noise_matched = noise_matched(correct_cells_SFTF,:,:);
noise_matched = noise_OD(iviv_vector(:,1)==1,:,:);

% get the cell IDs
cell_id = cat(1,str.cellID);
% get ids from the iviv vector
iviv_id = iviv_vector(iviv_vector(:,1)==1,2);

% for all the matched cells
for cells = 1:size(noise_matched,1)
    % get the boolean for selection of the cell
    id_bool = cell_id==iviv_id(cells);
    % if the cell is not here, print the id and skip
    if sum(id_bool)==0
        fprintf(strjoin({'Cell absent:',num2str(iviv_id(cells)),'\r\n'},'_'))
        continue
    end
%     % take the correlation for the preferred direction
%     [~,~,bin] = histcounts(str(id_bool).DIRpref,-22.5:45:382.5);
%     % rectify the last bin
%     if bin == 9
%         bin = 1;
%     elseif bin == 0
%         str(id_bool).noise = NaN;
%         continue
%     end
    % check the preference of the cell and take the corresponding average
    if str(id_bool).contra == 1
%         str(id_bool).noise = squeeze(noise_matched(cells,bin,1));
        str(id_bool).noise = squeeze(mean(noise_matched(cells,:,1),2));
    elseif str(id_bool).ipsi == 1
%         str(id_bool).noise = squeeze(noise_matched(cells,bin,2));
        str(id_bool).noise = squeeze(mean(noise_matched(cells,:,2),2));
    else
        if str(id_bool).ODIpref > 0
%             str(id_bool).noise = squeeze(noise_matched(cells,bin,1));
            str(id_bool).noise = squeeze(mean(noise_matched(cells,:,1),2));
        else
%             str(id_bool).noise = squeeze(noise_matched(cells,bin,2));
            str(id_bool).noise = squeeze(mean(noise_matched(cells,:,2),2));
        end
    end
end
%% Turn all the empty fields into NaN

% get the fields
field_list = fields(str);
% get the number of fields
num_fields = length(field_list);
% for all the fields
for f = 1:num_fields
    % find a prototype field
    for cells = 1:cell_num
        if ~isempty(str(cells).(field_list{f}))
            template_size = size(str(cells).(field_list{f}));
            break
        end
    end
    % for all the cells
    for cells = 1:cell_num
        if isempty(str(cells).(field_list{f}))
            str(cells).(field_list{f}) = NaN(template_size);
        end
    end
end
%% Save the structure (keep moving this down as features are added)
save(fullfile(stage5b_allinfo_path,'str_iviv.mat'),'str')