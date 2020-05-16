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

% define a vector of distance to try
distance_vector = 10:20:500;
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
distance_values = zeros(distance_number, target_num,2);
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

% for all the distances
for distance = distance_number:-1:1
    
    fprintf(strjoin({'Current distance:',num2str(distance),...
        'out of',num2str(distance_number),'\r\n'},'_'))
    % define the distance at which to cut
    distance_threshold = distance_vector(distance);
       
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
            [~,soma_idx] = min(abs(X_single(1,:)-target_somas(cells,1)));
            % determine the extent of the sum (balance depending on soma)
            if soma_idx <= size(X_single,2)/2
                extent = soma_idx-1;
            else
                extent = size(X_single,2) - soma_idx;
            end
            
            % isolate layer
            layer = interp_map(layer_idx1*grid_spacing:layer_idx2*grid_spacing,:);
            % get the horizontal input ratio
            fraction_values(cells) = sum(sum(layer(:,soma_idx-extent:soma_idx)))./...
                sum(sum(layer(:,soma_idx-extent:soma_idx+extent)));
        end
        
%         % determine the cutting position
%         cutting_location = target_somas(cells,1) + distance_threshold;
%         % find the position in the actual map
%         [~,cutting_map] = min(abs(X_single(1,:)-cutting_location));

        % get the center coordinates of the map

        center_x = target_somas(cells,1);
        center_y = -target_somas(cells,2);

        % get the logical mask
        mask = create_mask(X_single,Y_single,center_x,center_y,distance_threshold);
        
%         % plot mask and soma
%         [~,center_x_coord] = min(abs(X_single(1,:)-target_somas(cells,1)));
%         [~,center_y_coord] = min(abs(Y_single(:,1)-target_somas(cells,2)));
%         imagesc(mask)
%         hold on
%         plot(center_x_coord,size(Y_single,1)-center_y_coord,'ro')
        
        % get the map and layers
        cut_map = interp_map(layer_idx1*grid_spacing:layer_idx2*grid_spacing,:);
        % filter the map
        cut_map = cut_map.*mask(layer_idx1*grid_spacing:layer_idx2*grid_spacing,:);
   
%         % cut the map
%         cut_map(:,cutting_map:end) = 0;
%         cutting_location2 = target_somas(cells,1) - distance_threshold;
%         % find the position in the actual map
%         [~,cutting_map2] = min(abs(X_single(1,:) - cutting_location2));
%         % cut the map symmetrically
%         cut_map(:,1:cutting_map2) = 0;      

        % get the components
        cc = bwconncomp(cut_map);
        rp = regionprops(cc,cut_map,{'Area','WeightedCentroid'});
        % if there's more than one, leave the largest
        if cc.NumObjects > 1
            areas = cat(1,rp.Area);
            [~,idx] = max(areas);
            rp = rp(idx);
        elseif isempty(cc) || isempty(rp)
            continue
        end
        % store the centroid
        centroid_coord = rp.WeightedCentroid;
        cx = (X_single(1,round(centroid_coord(1))) - target_somas(cells,1))./grid_spacing;
        cy = (-Y_single(round(centroid_coord(2)),1) - row_shift*grid_spacing - ...
            target_somas(cells,2))./grid_spacing;
        % calculate the length and store
        centroid_vector(cells) = sqrt(cx^2+cy^2);
        % calculate the angle
        centroid_angle(cells) = rad2deg(atan2(cy,cx));
    end
    
    % if it's the full distance, establish as the centroid vector
    if distance == distance_number
        original_centroids = centroid_vector;
    end
    
%     % calculate and store the correlation and pvalue
%     [distance_correlation(distance,1),distance_correlation(distance,2)] =...
%         corr(centroid_vector,target_ori);
    % store the entire curve
    distance_values(distance,:,1) = centroid_vector;
    distance_values(distance,:,2) = centroid_angle;

    % for all the cell types
    for celltype = 1:celltype_num
        % get the idx vector
        idx = celltype_matrix(:,celltype);
        % calculate the correlation
%         [distance_percelltype(distance,celltype,1),distance_percelltype(distance,celltype,2)] = ...
%             corr(centroid_vector(idx),target_ori(idx));
        [distance_percelltype(distance,celltype,1),~] = ...
            corr(centroid_vector(idx),original_centroids(idx));
        distance_percelltype(distance,celltype,2) = mean(centroid_vector(idx));
    end
    
    
    
end
%% Plots per cell type

close all

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 300, 250]);
% mean_celltype = mean(distance_percelltype,1);
% std_celltype = std(distance_percelltype,0,1);
colors = {'k','b','r'};

% for all celltypes
for celltype = 1:celltype_num
    plot(distance_vector,distance_percelltype(:,celltype,1),colors{celltype})
%     plot(distance_percelltype(:,celltype,2),distance_percelltype(:,celltype,1),colors{celltype})
%     plot3(distance_vector,distance_percelltype(:,celltype,2),distance_percelltype(:,celltype,1),colors{celltype})
    hold on
scatter(distance_vector,distance_percelltype(:,celltype,1),distance_percelltype(:,celltype,2).*50,colors{celltype},'filled')
box off;    
%     plot(distance_vector,distance_percelltype(:,celltype,1),colors{celltype})
%     plot(distance_vector,distance_percelltype(:,celltype,2),strcat('--',colors{celltype}))
%     plot(distance_vector([1 distance_number]),[0.05 0.05],'m--')
end
% xlabel('Average vector length')
xlabel('Cutting point')
ylabel('Correlation')
set(gca,'FontSize',10)
%legend({'NR','Aligned','Ortho'},'location','southeast')
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


% plot the average +/- sem of the group centroid trajectories
figure
for celltype = 1:celltype_num
    % get the corresponding celltype set
    current_celltype = distance_values(:,celltype_matrix(:,celltype),:);
    % get the mean and sem
    mean_type = mean(current_celltype(:,:,1),2);
    sem_type = std(current_celltype(:,:,1),0,2)./sqrt(size(current_celltype,2));

    shadedErrorBar(distance_vector,mean_type,sem_type,'lineProps',colors{celltype})
    hold on
    
end

ylabel('Vector length')
xlabel('Cutting distance')