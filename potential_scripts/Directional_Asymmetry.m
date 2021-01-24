%% Work with the cells in the correct angle

%% Clean up

clearvars
close all
matlabrc
Paths
%% Load the relevant files

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Get the polarities

% allocate memory for the cells
polarity_cell = cell(2,1);
% for both setups
for setup = 1:2
    ori_setup = strcat('Setup',num2str(setup));
    polarity_cell{setup} = ...
        load(fullfile(preprocessing_path,strcat('ExcInh_map_',ori_setup,'.mat')),'ExcInh_map');
    polarity_cell{setup} = polarity_cell{setup}.ExcInh_map;
end

% concatenate them
polarity_vector = cat(1,polarity_cell{:});
%% Get the cells

angle_vector = cat(1,str.DIRpref);
selection_vector = cat(1,str.DSIpref)>0.25;

% define the target directions
target_directions = [130 310];
% target_directions = [45 225];

% define the width of the cone (in degrees)
cone_width = 30;

% create the vectors

% allocate memory for the vectors
target_vector = zeros(size(angle_vector,1),length(target_directions));

% for all the target directions
for dirs = 1:length(target_directions)
    % offset the directions to allow rotation
    temp_dirs = target_directions(dirs);% + 360;
    temp_angles = angle_vector;% + 360;
    % get the angles
    target_vector(:,dirs) = (temp_angles>(temp_dirs-cone_width))&...
        (temp_angles<(temp_dirs+cone_width))&selection_vector;
end

% concatenate
target_combined = any(target_vector,2); 

% get the cells
target_cells = str(target_combined);
% also get a vector for the direction side
direction_side = target_vector(target_combined,1);

% just leave 1 cell
target_cells = target_cells(1);
direction_side = direction_side(1);
% get the number of target cells
num_cells = length(target_cells);
%% Plot the input maps

close all

exc_fig = figure;
inh_fig = figure;
% for all the cells
for cells = 1:num_cells
    figure(exc_fig)
    subplot(round(sqrt(num_cells)),ceil(sqrt(num_cells)),cells)
    
    % get the delta map
%     delta_map = target_cells(cells).subpixel_excMap - target_cells(cells).subpixel_inhMap;
    delta_map = target_cells(cells).subpixel_raw_excMap;
    % flip the map if the it's the other direction
    if direction_side(cells)==0
        delta_map = fliplr(delta_map);
    end
    imagesc(delta_map)
     hold on
    plot([8.5 8.5],get(gca,'YLim'),'k')
    axis square
    set(gca,'TickLength',[0 0],'Xtick',[],'YTick',[])
    
    figure(inh_fig)
    subplot(round(sqrt(num_cells)),ceil(sqrt(num_cells)),cells)
    
    % get the delta map
%     delta_map = target_cells(cells).subpixel_excMap - target_cells(cells).subpixel_inhMap;
    delta_map = target_cells(cells).subpixel_raw_inhMap;
    % flip the map if the it's the other direction
    if direction_side(cells)==0
        delta_map = fliplr(delta_map);
    end
    imagesc(delta_map)
    hold on
    plot([8.5 8.5],get(gca,'YLim'),'k')

    axis square
    set(gca,'TickLength',[0 0],'Xtick',[],'YTick',[])
    
end
%% Plot the average delta map
close all
% allocate memory for the average maps
average_maps = zeros(16,16,2);
% for all the cells
for cells = 1:num_cells
    % get the maps
    excmap = target_cells(cells).subpixel_excMap;
    inhmap = target_cells(cells).subpixel_inhMap;
    % flip if needed
    if direction_side(cells) == 0
        excmap = fliplr(excmap);
        inhmap = fliplr(inhmap);
    end
    
    % add to the average
    average_maps(:,:,1) = average_maps(:,:,1) + excmap/num_cells;
    average_maps(:,:,2) = average_maps(:,:,2) + inhmap/num_cells;
end

% plot the maps
figure
subplot(2,2,1)
imagesc(average_maps(:,:,1))
hold on
plot([8.5 8.5],get(gca,'YLim'),'k')
subplot(2,2,3)
imagesc(average_maps(:,:,2))
hold on
plot([8.5 8.5],get(gca,'YLim'),'k')
% plot the average delta map
subplot(2,2,2)
imagesc(diff(average_maps,1,3))
hold on
plot([8.5 8.5],get(gca,'YLim'),'k')
%% Get the cell names

target_names = {target_cells.cellName};

target_ids = {target_cells.cellID};
%% Get the paths to all the maps

% define the target folders
folder_cell = {'Setup1','Setup2','Setup1_TTX','Setup2_TTX'};
path_info = get_map_paths(input_maps_path,folder_cell);
% concatenate
path_info = cat(1,path_info{:,1});


%get the number of maps to do
map_num = size(path_info,1);
%% Get an array with the list of cell names for each map

%extract the cell names from the folder all array
%allocate memory to save the cell names
cell_names = cell(map_num,1);
%for all the folders
for maps = 1:map_num
    temp_name = strsplit(path_info{maps},'\');
    cell_names{maps} = strcat(temp_name{end-2}(1:2),temp_name{end-3}(1:6),temp_name{end-2}(3:6));
end
%% Define a set of target names for testing

% target_names = {'SW1506160005','SW1606160002'};
%% Extract the target maps
tic
% get the indexes
% [~,ia] = intersect(cell_names,target_names);
ia = ismember(cell_names,target_names);

% get the subset of cell names
name_subset = cell_names(ia);

% get the paths
target_paths = path_info(ia);

% % get the polarities (cut the vector artificially until we figure out)
file_polarities = [polarity_vector{ia,1}]';

% % leave only 1 map
% name_subset = name_subset(1);
% target_paths = target_paths(1);
% file_polarities = file_polarities(1);

% get the number of files
num_files = length(target_paths);

% define the binning
binning = 10;

% allocate memory for the maps (arbitrary, if more flexible, need to setup
% varargin on the xsg_extractor)
maps_matrix = zeros(2000/binning,256,2,num_cells);
% allocate memory for the info and the soma for both polarities
info_cell = cell(num_cells,2,2);

% define the options
opts = struct([]);
opts(1).binning = binning;

% for all the indexes
for files = 1:num_files
    disp(files)
    % get the cell index
    cell_idx = find(contains(target_names,name_subset{files}));
    % get the total number of reps for this polarity and cell
    reps = sum((contains(name_subset,target_names{cell_idx}))&(file_polarities(files)==file_polarities));
    % get the map and info 
    % TODO: THIS ASSUMES STABILITY BETWEEN REPS, NEED TO IMPROVE THAT
    [current_map,info_cell{cell_idx,file_polarities(files)+1,1},...
        info_cell{cell_idx,file_polarities(files)+1,2}] = ...
        xsg_extractor(target_paths{files},file_polarities(files)+1,opts);
    
    % DEBUGGING
    % artificially mark a position as NaN
    if files < 3
        temp = info_cell{cell_idx,file_polarities(files)+1,1};
        temp(149,5) = 0;
        info_cell{cell_idx,file_polarities(files)+1,1} = temp;
    end
    
%     % interpolate the map
%     [current_map,info_cell{cell_idx,file_polarities(files)+1,1}] = ...
%         interpolate_map(info_cell{cell_idx,file_polarities(files)+1,1},current_map);
    
    % get the map
    maps_matrix(:,:,file_polarities(files)+1,cell_idx) = ...
        nansum(cat(5,maps_matrix(:,:,file_polarities(files)+1,cell_idx), ...
        current_map./reps),5);
end
toc
%% Calculate the actual maps

% calculate the overall centered soma for the subpixel interpolation
centered_soma = calculate_centered_soma(str);

% allocate memory for the maps
ready_maps = cell(num_cells,1);

% for all the cells
for cells = 1:num_cells
    
    % get the raw maps
    raw_maps = map_plotter(info_cell(cells,:,1),maps_matrix(:,:,:,cells));
    
    % get the onset maps
    onset_exc = reshape(info_cell{cells,1,1}(:,7),16,16);
    onset_inh = reshape(info_cell{cells,2,1}(:,7),16,16);
    % add the onset maps
    raw_maps = cat(3,raw_maps,onset_exc,onset_inh);
    
    % apply rotations, alignment and such
    
    % get the sliceOri
    % if 1, rotate
    sliceOri = target_cells(cells).sliceOri;
    if sliceOri == 0
        raw_maps = fliplr(raw_maps);
    end
    % rotate to align direction
    if direction_side(cells)==0
        raw_maps = fliplr(raw_maps);
    end
    
    % get the inputs for the alignment
    setup = target_cells(cells).cellID>47;
    soma_center = target_cells(cells).somaCenter;
    
    % align the map at the subpixel level
    raw_maps = align_subpixel(raw_maps,sliceOri,centered_soma,soma_center,setup);
    
    % normalize the maps
    % for all dims
    for dims = 1:5
        % get the map
        curr_map = raw_maps(:,:,dims);
        % define the normalization factor
        if dims == 1
            factor = min(curr_map(:));
        else
            factor = max(curr_map(:));
        end
        % normalize
        raw_maps(:,:,dims) = curr_map./factor;
    end
    
    
    
    % save the maps
    ready_maps{cells} = raw_maps;
    
end

% concatenate
ready_maps = cat(4,ready_maps{:});
% autoArrangeFigures
%% Plot individual maps for FS and non FS
close all

exc = figure;
slow = figure;
fast = figure;

for cells = 1:num_cells
    
    figure(exc)
    subplot(round(sqrt(num_cells)),ceil(sqrt(num_cells)),cells)
    imagesc(ready_maps(:,:,1,cells))
    axis equal
    axis tight
    
    figure(slow)
    subplot(round(sqrt(num_cells)),ceil(sqrt(num_cells)),cells)
    imagesc(ready_maps(:,:,2,cells))
    axis equal
    axis tight
    
    figure(fast)
    subplot(round(sqrt(num_cells)),ceil(sqrt(num_cells)),cells)
    imagesc(ready_maps(:,:,3,cells))
    axis equal
    axis tight
end
% autoArrangeFigures
%% Plot the average maps

close all

% for the 3 types of map
for maptype = 1:5
    figure
    imagesc(nanmean(ready_maps(:,:,maptype,:),4))
end
autoArrangeFigures
%% Plot the onset maps

close all

% for all the cells
for cells = 1:num_cells
    figure
    current_map = info_cell{cells,2,1};
    imagesc(reshape(current_map(:,7),16,16))
end
autoArrangeFigures
%% Visualize the maps
close all

% for all the cells
for cells = 1:num_cells
    figure
    % for both polarities
    for polarity = 1:2
        % initialize a position counter
        pos_count = 1;
        % get the current map and normalize
        curr_map = maps_matrix(:,:,polarity,cells);
        
        % get the color
        switch polarity
            case 1
                color = 'r';
                factor = -max(abs(curr_map(:)));
            case 2
                color = 'b';
                factor = -max(curr_map(:));
        end
        
        curr_map = curr_map./factor;
        % for all the positions
        for x_pos = 1:16
            for y_pos = 1:16
                
                % get the x and y position to plot
                x = linspace(0,0.8,size(maps_matrix,1))+x_pos;
                y = y_pos-1;
                
                % plot
                plot(x,curr_map(:,pos_count)+y,color)
                hold on
                plot([x(100),x(100)],[0 1]+y,'--k')
                % increment the counter
                pos_count = pos_count + 1;
            end
        end
    
    end
    set(gca,'TickLength',[0 0],'YDir','reverse')
    
end