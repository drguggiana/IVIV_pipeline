%% Calculate the max value per layer

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
celltype_colors = {'k','b','r'};

% get the maps
% selection_vector = vertcat(str.resp)==0;
selection_vector = ones(length(str),1)==1;
angle_span = 25;
% get vectors to identify all cell types
non_resp_cells = vertcat(str.resp)==0;
% threshold based on function
OS_cells = vertcat(str.OSIpref)>0.25;
% % get the TW
% tw = vertcat(str.Sigmapref);
% OS_cells = tw<prctile(tw,25);

% % select the type of map
% map_type = 'inh';


% define the cells groups
aligned_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+125)...
    & vertcat(str.ORIpref)>(125-angle_span));
ortho_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+45)...
    & vertcat(str.ORIpref)>(45-angle_span));

% generate a matrix with all the vectors
celltype_matrix = horzcat(non_resp_cells,aligned_cells,ortho_cells);
% get the number of cell types
celltype_num = size(celltype_matrix,2);


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
            row_shift = 2;
        case 2
            layer_idx1 = 6;
            layer_idx2 = 7;
            row_shift = 5;
        case 3
            layer_idx1 = 8;
            layer_idx2 = 10;
            row_shift = 7;
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
            
            
            % if the max is 1 (i.e. at the corner of the image), leave a NaN
%             if max_idx == 1 || max_y == 1 || max_y == size(cut_map,1)
%                 max_values(cells,maptype,layer,2:3) = NaN;
%                 continue
%             end
            % store the distances in both axis from the soma
            max_values(cells,maptype,layer,2) = (max_x - soma_x)/69;
            max_values(cells,maptype,layer,3) = (max_y - ...
                (size(Y_single,1)-soma_y-(layer_idx1-1)*grid_spacing))/69;
        end
    end
end
%% Plot the distributions of max values
close all
% figure
maptype_colors = {'r','b'};
for maptype = 1:2
    figure
    % set up a counter
    counter = 1;
    switch maptype
        case 1
            field1 = 'ex';
        case 2
            field1 = 'in';
    end
    for layer = 1:3
        switch layer
            case 1
                field2 = 'L23';
            case 2
                field2 = 'L4';
            case 3
                field2 = 'L5';
        end
        % get the maxima
        maxima = squeeze(max_values(:,maptype,layer,:));
        % get the corresponding centroids
        centroids = vertcat(str.(strcat('ang_',field1,field2)));
        
        % for the 2 coordinates
        for coord = 1:2
            subplot(3,2,counter)
            % get the centroid to soma
            centroid_coord = centroids(:,2+coord) - centroids(:,coord);
            
            % calculate the correlation
            [x,y] = nan_remover(maxima(:,coord+1),centroid_coord);
            [rho,pval] = corr(x,y);
            
            scatter(x,y,30,maptype_colors{maptype})
            hold on
            title(strjoin({field2,num2str(rho),num2str(pval)},' '),'Interpreter','None')
            axis square
            % update the counter
            counter = counter + 1;
        end
        
    end
end

autoArrangeFigures
%% Compare the centroid positions per celltype

close all
% figure
maptype_colors = {'r','b'};
for maptype = 1:2
    figure
    % set up a counter
    counter = 1;
    switch maptype
        case 1
            field1 = 'ex';
        case 2
            field1 = 'in';
    end
    for layer = 1:3
        switch layer
            case 1
                field2 = 'L23';
            case 2
                field2 = 'L4';
            case 3
                field2 = 'L5';
        end
        % get the maxima
        maxima = squeeze(max_values(:,maptype,layer,:));
        % get the corresponding centroids
%         centroids = vertcat(str.(strcat('ang_',field1,field2,'_nonweighted')));
        centroids = vertcat(str.(strcat('ang_',field1,field2)));

        
        % for all the celltypes
        for celltype = 1:3
        
            cell_idx = celltype_matrix(:,celltype);
%         % for the 2 coordinates
%         for coord = 1:2
            subplot(3,3,counter)
            % get the centroid to soma
            centroid_x = centroids(cell_idx,3) - centroids(cell_idx,1);
            centroid_y = centroids(cell_idx,4) - centroids(cell_idx,2);
            
            maxima_x = maxima(cell_idx,2);
            maxima_y = maxima(cell_idx,3);
            % calculate the correlation
            [x,y] = nan_remover(maxima_x,centroid_x);
            [rho,pval] = corr(x,y);
            
            scatter(maxima_x.*69,maxima_y.*69,30,celltype_colors{celltype})
            hold on
            scatter(centroid_x.*69,centroid_y.*69,30,celltype_colors{celltype},'filled')
            title(strjoin({field2,num2str(rho),num2str(pval)},' '),'Interpreter','None')
            set(gca,'YDir','reverse')
            axis equal
            % update the counter
            counter = counter + 1;
%         end
        end
        
    end
    suptitle(field1)
end

% autoArrangeFigures
%% 
close all

selection_vector = ones(147,1)==1;
% selection_vector = vertcat(str.OSIpref)>0.25;

figure
L23 = vertcat(str.ang_inL23);
L4 = vertcat(str.ang_inL4);
L23x = L23(selection_vector,3).*69;
L23y = L23(selection_vector,4).*69;

somax = L23(selection_vector,1).*69;
somay = L23(selection_vector,2).*69;

L4x = L4(selection_vector,3).*69;
L4y = L4(selection_vector,4).*69;

% subplot(1,2,1)
scatter(L23x,L23y)
% axis equal
% subplot(1,2,2)
hold on
scatter(L4x,L4y)
scatter(somax,somay)

axis equal
set(gca,'YDir','reverse')
% [L23x,somay] = nan_remover(vertcat(str.pialD),somay);
% [L23y,L4y] = nan_remover(L23y,L4y);
% [rho,pval] = corr(abs(L23x),somay)
% [rho,pval] = corr(L23y,L4y)


L23 = vertcat(str.ang_exL23);
L4 = vertcat(str.ang_exL4);
L23x = L23(selection_vector,3).*69;
L23y = L23(selection_vector,4).*69;

somax = L23(selection_vector,1).*69;
somay = L23(selection_vector,2).*69;

L4x = L4(selection_vector,3).*69;
L4y = L4(selection_vector,4).*69;

scatter(L23x,L23y)
scatter(L4x,L4y)

% for cells = 1:sum(selection_vector)
%     plot([somax(cells),L23x(cells)],[somay(cells),L23y(cells)],'k')
% end
% scatter(somax,somay)
axis equal
set(gca,'YDir','reverse')
% % MAXIMA
% L23 = horzcat(str.max_in)';
% L4 = horzcat(str.max_in)';
% L23x = L23(:,2);
% L23y = L23(:,3);
% 
% % somax = L23(:,1);
% % somay = L23(:,2);
% 
% L4x = L4(:,5);
% L4y = L4(:,6);
% 
% subplot(1,2,2)
% scatter(L23x,L23y)
% % axis equal
% % subplot(1,2,2)
% hold on
% scatter(L4x,L4y)
% scatter(somax,somay)
% 
% axis equal
% set(gca,'YDir','reverse')
% 
% 
% 
% L23 = horzcat(str.max_ex)';
% L4 = horzcat(str.max_ex)';
% L23x = L23(:,2);
% L23y = L23(:,3);
% 
% somax = L23(:,5);
% somay = L23(:,6);
% 
% L4x = L4(:,3);
% L4y = L4(:,4);
% 
% % subplot(1,2,2)
% scatter(L23x,L23y)
% % axis equal
% % subplot(1,2,2)
% hold on
% scatter(L4x,L4y)
% % scatter(somax,somay)
% axis equal
% set(gca,'YDir','reverse')

