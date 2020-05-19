function max_values = max_perlayer_calculation(str)
%% Attempt with cutting the cloud

% get the maps
% selection_vector = vertcat(str.resp)==0;
selection_vector = ones(length(str),1)==1;

% get the soma positions
target_somas = cat(1,str(selection_vector).subpixel_soma);

% define the grid spacing in microns
grid_spacing = 69;

% get the map size in microns
map_size = round(size(str(1).subpixel_excMap,1).*grid_spacing);
% get the map limits
map_lim = map_size/2-grid_spacing/2;

% create the grid
[Y,X] = ndgrid(-map_lim:grid_spacing:map_lim,...
    -map_lim:grid_spacing:map_lim);

% define the single micron grid
[Y_single,X_single] = ndgrid(-map_lim:map_lim,...
    -map_lim:map_lim);

% allocate memory for the max values and coordinates
max_values = zeros(sum(selection_vector),2,4,3);
% for all the distances
for layer = 1:4
    
    fprintf(strjoin({'Current layer:',num2str(layer),'\r\n'},' '))
    switch layer
        case 1
            layer_idx1 = 3;
            layer_idx2 = 5;
        case 2
            layer_idx1 = 6;
            layer_idx2 = 7;
        case 3
            layer_idx1 = 8;
            layer_idx2 = 10;
        case 4
            layer_idx1 = 1;
            layer_idx2 = 16;
    end
    
    % for both map types
    for maptype = 1:2
        
        switch maptype
            case 1
                map_type = 'exc';
            case 2
                map_type = 'inh';
        end
        % get the maps
        target_maps = cat(3,str(selection_vector).(strcat('subpixel_',map_type,'Map')));

        % get the number of maps
        target_num = size(target_maps,3);
        % for all the cells
        for cells = 1:target_num
            
            % get the corresponding map
            map = target_maps(:,:,cells);
            % interpolate the map
            interpolant = griddedInterpolant(Y,X,map);
            interp_map = interpolant(Y_single,X_single);
            
            % get the position of the soma in index
            [~,soma_x] = min(abs(X_single(1,:)-target_somas(cells,1)));
            [~,soma_y] = min(abs(Y_single(:,1)-target_somas(cells,2)));
            
            % get the layer boundaries
            [~,idx1] = min(abs(Y_single(:,1)-((layer_idx1-1)*grid_spacing-517.5)));
            [~,idx2] = min(abs(Y_single(:,1)-((layer_idx2-1)*grid_spacing-517.5)));
 
            % get the map and layers
            cut_map = interp_map(idx1:idx2,:);
            
            % calculate the position and value of the maximum and store
            [max_values(cells,maptype,layer,1),max_idx] = max(cut_map(:));
            % convert the coordinates into x and y
            [max_y,max_x] = ind2sub(size(cut_map),max_idx);
            
            % store the distances in both axis from the soma
            max_values(cells,maptype,layer,2) = (max_x - soma_x)/69;
            max_values(cells,maptype,layer,3) = (max_y - ...
                (size(Y_single,1)-soma_y-(layer_idx1-1)*grid_spacing))/69;
        end
    end
end