%% Clean up
clearvars
close all
Paths
% define the path to the structure
str_path = structure_file_path;
% load the structure
load(str_path)
%% Define the groups of cells

% define the plot colors (for celltypes)
colors = {'k','b','r'};
% define the names
morphos_name = {'Apical','Basal'};
layer_name = {'L23','L4'};
celltype_name = {'NR','Aligned','Ortho'};
% define the span of the angle cone
angle_span = 25;
% get vectors to identify all cell types
% non_resp_cells = vertcat(str.resp)==0&vertcat(str.sftf_resp)==0;
non_resp_cells = vertcat(str.resp)==0;
% threshold based on function
OS_cells = vertcat(str.OSIpref)>0.25;
% % get the TW
% tw = vertcat(str.Sigmapref);
% OS_cells = tw<prctile(tw,25);

% define the cells groups
aligned_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+135)...
    & vertcat(str.ORIpref)>(135-angle_span));
ortho_cells = OS_cells & (vertcat(str.ORIpref)<(angle_span+45)...
    & vertcat(str.ORIpref)>(45-angle_span));
% generate a matrix with all the vectors
celltype_matrix = horzcat(non_resp_cells,aligned_cells,ortho_cells);
% get the number of cell types
celltype_num = size(celltype_matrix,2);
%% Calculate the amounts of input that go with the morphology

% allocate memory for the results
input_morpho_overlap = cell(2,2);

% for exc and inh
for maptype = 1:2
    switch maptype
        case 1
            map_name = 'subpixel_excMap';
        case 2
            map_name = 'subpixel_inhMap';
    end
    % load the corresponding maps
    current_maps = cat(3,str.(map_name));
    % for apical and basal
    for morphotype = 1:2
        switch morphotype
            case 1
                morpho_name = 'morphoMap_apical_aligned';
            case 2
                morpho_name = 'morphoMap_basal_aligned';
        end
        % load the corresponding morpho
        current_morpho = cat(3,str.(morpho_name));
        % normalize morpho
        current_morpho = reshape(current_morpho,256,[]);
        % nan the 0s
        current_morpho(current_morpho==0) = NaN; 
        current_morpho = (current_morpho - nanmean(current_morpho,1))./nanstd(current_morpho,0,1);
        current_morpho = reshape(current_morpho,16,16,[]);
        % multiply the two
        combination = current_maps.*current_morpho;
        % turn 0 into NaN
        combination(combination==0) = NaN;
        % save the result
        input_morpho_overlap{maptype,morphotype} = combination;
        
    end
end
%% Plot the results
close all

% figure
% % initialize a counter
% counter = 1;
% % for exc and inh
% for maptype = 1:2
%     % for apical and basal
%     for morphotype = 1:2
%         subplot(2,2,counter)
%         
%         histogram(input_morpho_overlap{maptype,morphotype}(:),50)
%         
%         % update the counter
%         counter = counter + 1;
%         
%     end
% end

figure
% allocate memory for the averages
mean_combo = zeros(4,2,celltype_num);

% for all the groups
for celltype = 1:celltype_num
    % initialize a counter
    counter = 1;
    % for exc and inh
    for maptype = 1:2
        % for apical and basal
        for morphotype = 1:2
           
            % get the maps
            maps = input_morpho_overlap{maptype,morphotype}(:,:,celltype_matrix(:,celltype));

            mean_combo(counter,1,celltype) = nanmean(maps(:));
            mean_combo(counter,2,celltype) = nanstd(maps(:))./sqrt(size(maps,3));

            % update the counter
            counter = counter + 1;
            
        end
    end
    errorbar(1:counter-1,mean_combo(:,1,celltype),mean_combo(:,2,celltype),'o')
    hold on
end

set(gca,'XLim',[0 counter])
%% Calculate the morpho centroids

% allocate memory for the results
morpho_cell = cell(3,2);


% for apical and basal
for morphotype = 1:2
    switch morphotype
        case 1
            morpho_name = 'morphoMap_apical_aligned';
        case 2
            morpho_name = 'morphoMap_basal_aligned';
    end
    % load the corresponding morpho
    current_morpho = cat(3,str.(morpho_name));
    % get the soma positions
    soma_center = vertcat(str.subpixel_soma);
    
    % for all the layers
    for layer = 1:3
        switch layer
            case 1
                row_shift = 2;
                idx_start = 3;
                idx_end = 5;
            case 2
                row_shift = 5;
                idx_start = 6;
                idx_end = 7;
            case 3
                row_shift = 7;
                idx_start = 8;
                idx_end = 10;
        end
        % calculate the centroids
        out_ang = centroid_map(current_morpho(idx_start:idx_end,:,:),...
            soma_center(:,1),soma_center(:,2),cat(1,str.pialD),row_shift);
        % get the corrected centroids
        cx = abs(out_ang(:,3) - out_ang(:,1));
        cy = abs(out_ang(:,4) - out_ang(:,2));
        cvl = out_ang(:,8);
        cang = out_ang(:,5);
        
        %save the centroid position and vector
        morpho_cell{layer,morphotype} = horzcat(cx,cy,cvl,cang);
        
    end
    
end
%% Plot the centroid results

close all

figure
% initialize a plot counter
counter = 1;
% for all the layers
for layer = 1:2
    
    % for apical and basal
    for morphotype = 1:2
        
        subplot(2,2,counter)
        % plot the centroid
        current_data = morpho_cell{layer,morphotype};
        % for all the groups
        for celltype = 1:celltype_num
            scatter(current_data(celltype_matrix(:,celltype),1),...
                current_data(celltype_matrix(:,celltype),2),colors{celltype})
            hold on
            
        end
        title(strjoin({morphos_name{morphotype},layer_name{layer}},' '))
        % update the counter
        counter = counter + 1;
    end
end

figure

% initialize a plot counter
counter = 1;
% for all the layers
for layer = 1:2
    
    % for apical and basal
    for morphotype = 1:2
        
        subplot(2,2,counter)
        % plot the centroid
        current_data = morpho_cell{layer,morphotype};
        % allocate memory to store the averages
        temp_average = zeros(celltype_num,2);
        % for all the groups
        for celltype = 1:celltype_num
            temp_average(celltype,1) = nanmean(current_data(celltype_matrix(:,celltype),1));    
            temp_average(celltype,2) = nanstd(current_data(celltype_matrix(:,celltype),1))./sqrt();
        end
        bar(1:celltype_num,temp_average(:,1))
        hold on
        errorbar(1:celltype_num,temp_average(:,1),temp_average(:,2),'o')
        title(strjoin({morphos_name{morphotype},layer_name{layer}},' '))
        % update the counter
        counter = counter + 1;
    end
end

autoArrangeFigures
%% Relate to input
close all
main = figure;
pval = figure;

% initialize a plot counter
counter = 1;

% for all the layers
for layer = 1:2
    switch layer
        case 1
            field_name = 'ang_exL23';
        case 2
            field_name = 'ang_exL4';
        
    end
    % get the input maps
    input_cells = vertcat(str.(field_name));
%     input_cells = input_cells(:,8);
    input_cells = abs(input_cells(:,4)-input_cells(:,2));

    % for apical and basal
    for morphotype = 1:2
        
        % plot the centroid
        current_data = morpho_cell{layer,morphotype};
        % allocate memory to store the averages
        temp_corr = zeros(celltype_num,2);
        % for all the groups
        for celltype = 1:celltype_num
%             temp_average(celltype,1) = nanmean(current_data(celltype_matrix(:,celltype),3));    
%             temp_average(celltype,2) = nanstd(current_data(celltype_matrix(:,celltype),3));
           % get the morpho cells
           morpho_cells = current_data(celltype_matrix(:,celltype),2);
           % get the nanvector
           nan_vector = ~isnan(morpho_cells);
           % filter the morpho
           morpho_cells = morpho_cells(nan_vector);
           % get the input and filter
           input_temp = input_cells(celltype_matrix(:,celltype));
           input_temp = input_temp(nan_vector);
           
           if isempty(morpho_cells)
               continue
           end
           % calculate the correlation
           [temp_corr(celltype,1),temp_corr(celltype,2)] = corr(morpho_cells,input_temp);
           sum(nan_vector) 
        end
        figure(main)
        subplot(2,2,counter)
        
        scatter(1:celltype_num,temp_corr(:,1),30,[0 0 0;0 0 1;1 0 0],'o')
        hold on
        title(strjoin({morphos_name{morphotype},layer_name{layer}},' '))
        set(gca,'XLim',[0 celltype_num+1],'XTick',1:celltype_num,...
            'XTickLabels',celltype_name,'XTickLabelRotation',45)
        ylabel('Correlation')
        figure(pval)
        subplot(2,2,counter)
        
        scatter(1:celltype_num,temp_corr(:,2),30,[0 0 0;0 0 1;1 0 0],'o')
        set(gca,'YScale','log')
        hold on
        ylabel('p value')
        plot([0 celltype_num+1],[0.05 0.05],'m')
        %         errorbar(1:celltype_num,temp_average(:,1),temp_average(:,2),'o')
        title(strjoin({morphos_name{morphotype},layer_name{layer}},' '))
        set(gca,'XLim',[0 celltype_num+1],'XTick',1:celltype_num,...
            'XTickLabels',celltype_name,'XTickLabelRotation',45)
        % update the counter
        counter = counter + 1;
    end
end
%% Plot a property over the centroids
close all

% define the color code
color_code = vertcat(str.ORIpref);

figure
% initialize a plot counter
counter = 1;
% for all the layers
for layer = 1:2
    
    % for apical and basal
    for morphotype = 1:2
        
        subplot(2,2,counter)
        % plot the centroid
        current_data = morpho_cell{layer,morphotype};

        scatter(current_data(:,1),...
            current_data(:,2),30,color_code)
        hold on
            
        title(strjoin({morphos_name{morphotype},layer_name{layer}},' '))
        % update the counter
        counter = counter + 1;
    end
    axis equal
end