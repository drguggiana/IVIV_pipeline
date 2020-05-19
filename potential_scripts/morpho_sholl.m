%% Calculate the input map maxima appearing in concentric rings around the soma
%% Clean up
clearvars
close all
Paths
% define the path to the structure
str_path = structure_file_path;
% load the structure
load(str_path)
%% Attempt with cutting the cloud
close all
% define the celltype colors
colors = {'k','b','r'};

% get the maps
%     selection_vector = vertcat(str.OSIpref)>0.25;
%     selection_vector = selection_vector & (vertcat(str.ORIpref)<(angle_span+target_angle)...
%         & vertcat(str.ORIpref)>(angle_span-target_angle));
% selection_vector = vertcat(str.resp)==0;
selection_vector = ones(length(str),1)==1;
angle_span = 25;
% get vectors to identify all cell types
% non_resp_cells = vertcat(str.resp)==0&vertcat(str.sftf_resp)==0;
non_resp_cells = vertcat(str.resp)==0;
% threshold based on function
OS_cells = vertcat(str.OSIpref)>0.25;
% % get the TW
% tw = vertcat(str.Sigmapref);
% OS_cells = tw<prctile(tw,25);

% select the type of map
map_type = 'inh';

% select the layer
layer = 23;
switch layer
    case 23
        layer_idx1 = 3;
        layer_idx2 = 5;
        row_shift = 2;
    case 4
        layer_idx1 = 6;
        layer_idx2 = 7;
        row_shift = 5;
    case 5
        layer_idx1 = 8;
        layer_idx2 = 10;
        row_shift = 7;
    case 'all'
        layer_idx1 = 1;
        layer_idx2 = 16;
end
% define the cells groups
aligned_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+125)...
    & vertcat(str.ORIpref)>(125-angle_span));
ortho_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+45)...
    & vertcat(str.ORIpref)>(45-angle_span));

% generate a matrix with all the vectors
celltype_matrix = horzcat(non_resp_cells,aligned_cells,ortho_cells);
% get the number of cell types
celltype_num = size(celltype_matrix,2);

% get the maps
target_maps = cat(3,str(selection_vector).(strcat('subpixel_',map_type,'Map')));

% define the ring effective diameter
ring_diameter = 5;
% define a vector of distance to try
distance_vector = 10:ring_diameter:125;
distance_vector = horzcat(distance_vector,8*69);
% get the number of distances
distance_number = length(distance_vector);
% % allocate memory to store the results
% distance_correlation = zeros(distance_number,2);
% allocate memory to calculate the correlation per group
distance_percelltype = zeros(distance_number,celltype_num,2);

% get the number of maps
target_num = size(target_maps,3);

% allocate memory to store the actual asymmetries
distance_values = zeros(distance_number, target_num,3);
% allocate memory for the horizontal fraction of input
fraction_values = zeros(target_num,1);

% get the soma positions
target_somas = cat(1,str(selection_vector).subpixel_soma);
% get the orientations
target_ori = vertcat(str(selection_vector).ORIpref);

% get the number of maps
number_of_targets = size(target_maps,3);

% define the grid spacing in microns
grid_spacing = 69;

% get the map size in microns
map_size = round(size(target_maps,1).*grid_spacing);
% get the map limits
map_lim = map_size/2-grid_spacing/2;

% create the grid
[Y,X] = ndgrid(-map_lim:grid_spacing:map_lim,...
    -map_lim:grid_spacing:map_lim);

% define the single micron grid
[Y_single,X_single] = ndgrid(-map_lim:map_lim,...
    -map_lim:map_lim);

% allocate memory to store the masks
mask_cell = cell(target_num,2);

% initialize a frame counter
frame_counter = 1;
% define the path for the gif

% for all the distances
for distance = distance_number:-1:2
    
    fprintf(strjoin({'Current distance:',num2str(distance),...
        'out of',num2str(distance_number),'\r\n'},'_'))
    % define the distance at which to cut
    distance_threshold = distance_vector(distance-1);
       
    % get the slice ori cloud centers for setup 2
    centroid_vector = zeros(target_num,1);
    centroid_angle = zeros(target_num,1);
    
    % for all the cells
    for cells = 1:target_num
        
        % get the corresponding map
        map = target_maps(:,:,cells);
        % interpolate the map
        interpolant = griddedInterpolant(Y,X,map);
        interp_map = interpolant(Y_single,X_single);
        
        % if it's the full map, calculate the input fraction horizontally
        if distance == distance_number
            % get the position of the soma in index
            [~,soma_x] = min(abs(X_single(1,:)-target_somas(cells,1)));
            [~,soma_y] = min(abs(Y_single(:,1)-target_somas(cells,2)));
            % determine the extent of the sum (balance depending on soma)
            if soma_x <= size(X_single,2)/2
                extent = soma_x-1;
            else
                extent = size(X_single,2) - soma_x;
            end
            % get the layer boundaries
            [~,idx1] = min(abs(Y_single(:,1)-((layer_idx1-1)*grid_spacing-517.5)));
            [~,idx2] = min(abs(Y_single(:,1)-((layer_idx2-1)*grid_spacing-517.5)));
            
            % isolate layer
%             layer = interp_map(layer_idx1*grid_spacing:layer_idx2*grid_spacing,:);
            layer = interp_map(idx1:idx2,:);

            % get the horizontal input ratio
            fraction_values(cells) = sum(sum(layer(:,soma_x-extent:soma_x)))./...
                sum(sum(layer(:,soma_x-extent:soma_x+extent)));
        end

        % get the center coordinates of the map
        center_x = target_somas(cells,1);
        center_y = -target_somas(cells,2);
        
        % if it's the first run, calculate an additional mask
        if distance == distance_number
            mask_cell{cells,1} = create_mask(X_single,Y_single,...
                center_x,center_y,distance_vector(distance_number));
        else
            % otherwise overwrite the old one
            mask_cell{cells,1} = mask_cell{cells,2};
        end

        % get the logical mask
        mask_cell{cells,2} = create_mask(X_single,Y_single,center_x,center_y,distance_threshold);
        
        % produce the annular mask
        mask = xor(mask_cell{cells,1},mask_cell{cells,2});

        % get the map and layers
        cut_map = interp_map(idx1:idx2,:);
        % filter the map
        cut_map = cut_map.*mask(idx1:idx2,:);
        % calculate the position and value of the maximum and store
        [distance_values(distance,cells,1),max_idx] = max(cut_map(:));
        % convert the coordinates into x and y
        [max_y,max_x] = ind2sub(size(cut_map),max_idx);
        
%         % plot mask and soma
%         [~,center_x_coord] = min(abs(X_single(1,:)-target_somas(cells,1)));
%         [~,center_y_coord] = min(abs(Y_single(:,1)-target_somas(cells,2)));
%         imagesc(mask)
%         hold on
%         plot(center_x_coord,size(Y_single,1)-center_y_coord,'ro')

%         % make a gif of the rings
%         imagesc(cut_map)
%         hold on
%         plot(max_x,max_y,'ro')
%         axis equal
%         pause(0.1)
%         make_gif_frame(gcf,frame_counter,fullfile(gif_path,'test.gif'),0.1)
%         frame_counter = frame_counter + 1;
        
        % if the max is 1 (i.e. at the corner of the image), leave a NaN
        if max_idx == 1 || max_y == 1 || max_y == size(cut_map,1)
            distance_values(distance,cells,2:3) = NaN;
            continue
        end
        % store the distances in both axis from the soma
        distance_values(distance,cells,2) = max_x - soma_x;
        distance_values(distance,cells,3) = max_y - ...
            (size(Y_single,1)-soma_y-(layer_idx1-1)*grid_spacing);
    end
end

% remove the first row with only 0s
distance_values = distance_values(2:end,:,:);
distance_vector = distance_vector(2:end);
%% Plot the distributions of max values
close all
% figure
% define the offset
offset = 500;
% for all the cell types
for celltype = 1:celltype_num
%     subplot(1,3,celltype)
    figure
    % get the idx vector
    idx = celltype_matrix(:,celltype);
    % plot value and location of the maxima
    current_val = distance_values(:,idx,:);
    % get a colormap for the cells
    cmap = lines(size(current_val,2));
    % for all the cells
    for cells = 1:size(current_val,2)
        
        current_cell = current_val(:,cells,:);
        %     histogram(current_val(20,:,3),30)
        %     current_val = reshape(current_val,[],3);
        % %     current_val = mean(current_val,2);
        %     current_val = squeeze(current_val(10,:,:));
        scatter(current_cell(:,2)+offset*(cells-1),current_cell(:,3),30,current_cell(:,1))
        hold on
        plot(current_cell(:,2)+offset*(cells-1),current_cell(:,3),'color',cmap(cells,:))
    end
    set(gca,'YDir','reverse')
    axis equal
    
end
autoArrangeFigures
%% Plot the maxima linearly

close all
% for all the cell types
for celltype = 1:celltype_num
%     subplot(1,3,celltype)
    figure
    % get the idx vector
    idx = celltype_matrix(:,celltype);
    % plot value and location of the maxima
    current_val = distance_values(:,idx,:);
    
    % plot the values as a matrix
    subplot(1,2,1)
    imagesc(current_val(:,:,1))
    
    % plot the euclidean displacement
    euc_displacement = sqrt(diff(current_val(:,:,2),1,1).^2 + diff(current_val(:,:,3),1,1).^2);
    subplot(1,2,2)
    imagesc(log(euc_displacement))
    
%     % for all the cells
%     for cells = 1:size(current_val,2)
%         
%         current_cell = current_val(:,cells,:);
%         %     histogram(current_val(20,:,3),30)
%         %     current_val = reshape(current_val,[],3);
%         % %     current_val = mean(current_val,2);
%         %     current_val = squeeze(current_val(10,:,:));
%             
% 
%         
% %         scatter(current_cell(:,2)+offset*(cells-1),current_cell(:,3),30,current_cell(:,1))
% %         hold on
% %         plot(current_cell(:,2)+offset*(cells-1),current_cell(:,3),'color',cmap(cells,:))
%     end
%     set(gca,'YDir','reverse')
%     axis equal
    
end
autoArrangeFigures
%% Plots per cell type

close all

figure
% mean_celltype = mean(distance_percelltype,1);
% std_celltype = std(distance_percelltype,0,1);

% for all celltypes
for celltype = 1:celltype_num
    plot(distance_vector,distance_percelltype(:,celltype,1),colors{celltype})
%     plot(distance_percelltype(:,celltype,2),distance_percelltype(:,celltype,1),colors{celltype})
%     plot3(distance_vector,distance_percelltype(:,celltype,2),distance_percelltype(:,celltype,1),colors{celltype})
%     scatter(distance_vector,distance_percelltype(:,celltype,1),distance_percelltype(:,celltype,2).*50,colors{celltype})
    hold on
%     plot(distance_vector,distance_percelltype(:,celltype,1),colors{celltype})
%     plot(distance_vector,distance_percelltype(:,celltype,2),strcat('--',colors{celltype}))
%     plot(distance_vector([1 distance_number]),[0.05 0.05],'m--')
end
% xlabel('Average vector length')
xlabel('Cutting point')
ylabel('Correlation')
legend({'NR','Aligned','Ortho'},'location','southeast')
%% Plot fractions
close all
figure

% allocate memory for the bar values
fraction_mean = zeros(celltype_num,2);
% extract the averages per cell type
for celltype = 1:celltype_num
    
    fraction_mean(celltype,1) = mean(fraction_values(celltype_matrix(:,celltype)));
    fraction_mean(celltype,2) = std(fraction_values(celltype_matrix(:,celltype)))./...
        sqrt(sum(celltype_matrix(:,celltype)));
end

bar(fraction_mean(:,1))
hold on
errorbar(1:celltype_num,fraction_mean(:,1),fraction_mean(:,2),'o')

ylabel('Horizontal layer 2/3 input fraction')
set(gca,'XTick',1:celltype_num,'XTickLabels',{'NR','Aligned','Ortho'})
%% Plots using all cells
close all
% plot all of the centroid trajectories for the aligned, ortho and nr groups
figure

for celltype = 1:celltype_num
    % get the current cells
    current_celltype = distance_values(:,celltype_matrix(:,celltype),:);
    
%     current_celltype(:,current_celltype(1,:,2)<0,2) = ...
%         current_celltype(:,current_celltype(1,:,2)<0,2) + 180;
%     current_celltype(:,current_celltype(end,:,2)<90,2) = ...
%         180 - current_celltype(:,current_celltype(end,:,2)<90,2);
    
    subplot(1,2,1)
    plot(distance_vector,current_celltype(:,:,1),colors{celltype})
    hold on
    subplot(1,2,2)
    plot(distance_vector,current_celltype(:,:,2),strcat(colors{celltype},'--'))
    hold on

end
subplot(1,2,1)
ylabel('Vector length')
xlabel('Cutting distance')

subplot(1,2,2)
ylabel('Vector angle')
xlabel('Cutting distance')
% plot the curves in an image
figure
[~,sort_idx] = sort(vertcat(str.pialD));
imagesc(distance_values(:,sort_idx,1)')
xlabel('Cutting distance')
ylabel('Cells')
cba = colorbar;
ylabel(cba,'Vector length')