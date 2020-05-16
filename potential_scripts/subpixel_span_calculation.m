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
colors = {'k','b','r'};

% get vectors to identify all cell types
% non_resp_cells = vertcat(str.resp)==0&vertcat(str.sftf_resp)==0;
non_resp_cells = vertcat(str.resp)==0;
% threshold based on function
OS_cells = vertcat(str.OSIpref)>0.25;
% % get the TW
% tw = vertcat(str.Sigmapref);

% allocate memory for the horizontal fraction of input
subpixel_span = zeros(size(OS_cells,1),2);
% OS_cells = tw<prctile(tw,25);
for maptype = 1:2
    switch maptype
        case 1
            map_type = 'exc';
        case 2
            % select the type of map
            map_type = 'inh';
    end

    % select the layer
    layer = 5;
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

    % get the number of maps
    target_num = size(target_maps,3);

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

    % for all the cells
    for cells = 1:target_num

        % get the corresponding map
        map = target_maps(:,:,cells);
        % interpolate the map
        interpolant = griddedInterpolant(Y,X,map);
        interp_map = interpolant(Y_single,X_single);

        % get the x position of the soma in index
        [~,soma_idx] = min(abs(X_single(1,:)-target_somas(cells,1)));

        % isolate layer
        layer = interp_map(layer_idx1*grid_spacing:layer_idx2*grid_spacing,:);
        

        % get the span of the cloud
        % allocate memory for the temp span
        temp_span = zeros(2,1);
        % for both halves of the map
        for half = 1:2
            % cut the map at the soma
            switch half
                case 1
                    % left half
                    half_map = layer(:,1:soma_idx);
                    direction = 'last';
                    
                    % get the distance between the soma x position and the edge of the map
                    span_thres = max(prctile(half_map,10,2));
                    max_idx = soma_idx-find(max(half_map,[],1)<=span_thres,1,direction);
                case 2
                    % right half
                    half_map = layer(:,soma_idx:end);
                    direction = 'first';
                    % get the distance between the soma x position and the edge of the map
                    span_thres = max(prctile(half_map,10,2));
                    max_idx = find(max(half_map,[],1)<=span_thres,1,direction);
            end

            if isempty(max_idx)
                continue
            end
            temp_span(half) = max_idx;
        end
        % save the max distance
        subpixel_span(cells,maptype) = sum(temp_span);

    end
end
%% Plot the relationship with orientation

close all

for maptype = 1:2
    figure
    % for all celltypes
    for celltype = 1:celltype_num
    %     scatter(target_param(celltype_matrix(:,celltype)),subpixel_span(celltype_matrix(:,celltype)))

    %     [N,edges] = histcounts(subpixel_span(celltype_matrix(:,celltype)),'Normalization','probability');
    %     plot(N,colors{celltype})
        errorbar(celltype,mean(subpixel_span(celltype_matrix(:,celltype),maptype)),...
            std(subpixel_span(celltype_matrix(:,celltype),maptype))./sqrt(sum(celltype_matrix(:,celltype))),...
            colors{celltype},'marker','o')
        hold on

    end
    set(gca,'XLim',[0 celltype+1])
    set(gca,'XTick',1:celltype_num,'XTickLabels',{'NA','Aligned','Ortho'})
end

figure
target_param = vertcat(str.ORIpref);
target_param = target_param(vertcat(str.OSIpref)>0.25);
sub_span = subpixel_span(vertcat(str.OSIpref)>0.25,:);

scatter(sub_span(:,1),target_param)
[cor,pval] = circ_corrcl(deg2rad(target_param),sub_span(:,1));
hold on
scatter(sub_span(:,2),target_param)
[cor2,pval2] = circ_corrcl(deg2rad(target_param),sub_span(:,2));
title(strjoin({num2str(cor),num2str(pval),num2str(cor2),num2str(pval2)},'_'),'Interpreter','None')