%% Simulate the position of input maps and get the possible centroid positions
%% Clean up
clearvars
close all
Paths
% define the path to the structure
str_path = structure_file_path;
% load the structure
load(str_path)
%% Generate new distributions of centroids by arbitrarily translating the maps

% define the number of shuffles
shuffle_number = 100;
% define the number of layers
layer_num = 3;
% get the cell number
cell_num = size(str,1);
% allocate memory to store the shuffled information (shuff, cell, exin, layer,
% xy)
shuffle_results = zeros(shuffle_number,cell_num,2,layer_num,2);
pialD = vertcat(str.pialD);

% get the somas
soma_coord = cat(1,str.subpixel_soma);

% for all the shuffles
for shuff = 1:shuffle_number
    fprintf(strjoin({'Current shuffle:',num2str(shuff),...
        'out of',num2str(shuffle_number),'\r\n'},' '))

    % for exin
    for maptype = 1:2
        
        switch maptype
            case 1
                % get the maps
                maps = cat(3,str.subpixel_excMap);
            case 2
                % get the maps
                maps = cat(3,str.subpixel_inhMap);
        end
        
        % get a set of random vertical displacements
        displacement = randi(6,cell_num,1)-3;
%         displacement = zeros(cell_num,1);

        % displace the maps
        % allocate memory for the maps
        displaced_maps = zeros(size(maps));
        % for all the cells
        for cells = 1:cell_num
            % displace the map
            displaced_maps(:,:,cells) = circshift(maps(:,:,cells),displacement(cells),1);
            % eliminate the portions of the map that rolled over
            displaced_maps(1:displacement(cells),:,cells) = 0;
        end
        % for all the layers
        for layer = 1:layer_num
            switch layer
                case 1
                    layer_idx1 = 2;
                    layer_idx2 = 6;
                    row_shift = 2;
                case 2
                    layer_idx1 = 5;
                    layer_idx2 = 8;
                    row_shift = 5;
                case 3
                    layer_idx1 = 7;
                    layer_idx2 = 11;
                    row_shift = 7;
                case 4
                    layer_idx1 = 1;
                    layer_idx2 = 16;
                    row_shift = 0;
            end
            % cut the maps
            cut_maps = displaced_maps(layer_idx1:layer_idx2,:,:);
            % get the centroid and store
            ang_out = centroid_map(cut_maps,soma_coord(:,1),soma_coord(:,2),pialD,row_shift);
            % store in the matrix
            shuffle_results(shuff,:,maptype,layer,:) = ang_out(:,[3,4]);
        end
    end

end
%% Plot the results

close all



exin_labels = {'Exc','Inh'};
% for the maptype
for maptype = 1:2
    switch maptype
        case 1
            field1 = 'ex';
        case 2
            field1 = 'in';
    end
    dots = figure;
    histo = figure;
    % load the real data
    % allocate memory for the data (cell, layer, x y)
    real_data = zeros(cell_num,3,2);
    % for all the layers
    for layer = 1:3
        switch layer
            case 1
                field2 = 'L23';
            case 2
                field2 = 'L4';
            case 3
                field2 = 'L5';
        end
        % get the data
        temp_data = vertcat(str.(strcat('ang_',field1,field2)));
        % load in the matrix
        real_data(:,layer,:) = temp_data(:,[3 4]);
    end
    
 
    % convert to microns
    temp_data = real_data.*69;
    % for all layers
    for layer = 1:layer_num
        
        switch layer
            case 1
                field2 = 'L23';
            case 2
                field2 = 'L4';
            case 3
                field2 = 'L5';
        end
        % get the coordinates and flatten
        x_coord = reshape(shuffle_results(:,:,maptype,layer,1),[],1); 
        y_coord = reshape(shuffle_results(:,:,maptype,layer,2),[],1);
        

    %     % convert to microns
        x_coord = x_coord.*69;
        y_coord = y_coord.*69;
  
        
        % plot
        figure(dots)
        scatter(x_coord,y_coord)
        hold on
        scatter(temp_data(:,layer,1),temp_data(:,layer,2))
        
        axis equal
        set(gca,'YDir','reverse')
        title(exin_labels{maptype})
        
        figure(histo)
        subplot(1,3,layer)
        histogram(temp_data(:,layer,2),'Normalization','pdf')
        hold on
        histogram(y_coord,'Normalization','pdf')
        
        [flag,pval] = kstest2(temp_data(:,layer,2),y_coord);
        title(strjoin({field1,field2,num2str(pval)},' '))
        legend({'Real Cy','Shuffle'})
    end
%     scatter(somas(:,1),piald)

end
%% Plot the soma and data
% close all
figure
ori_pref = vertcat(str.ORIpref);
% for the maptype
for maptype = 1:2
    switch maptype
        case 1
            field1 = 'ex';
        case 2
            field1 = 'in';
    end
    
%     figure
    % for all the layers
    for layer = 1
        switch layer
            case 1
                field2 = 'L23';
            case 2
                field2 = 'L4';
            case 3
                field2 = 'L5';
        end
        % get the data
        temp_data = vertcat(str.(strcat('ang_',field1,field2)));
        % load in the matrix
%         real_data(:,layer,:) = temp_data(:,[3 4]);
        
        % get the somas
        if layer == 1
            somas = temp_data(:,[1 2]);
            piald = vertcat(str.pialD)./69;
        end
        
%         scatter(temp_data(:,3)-temp_data(:,1),temp_data(:,4)-temp_data(:,2),30,ori_pref,'filled')
        scatter(temp_data(:,2),temp_data(:,4),30,ori_pref,'filled')
        xlabel('soma')
        ylabel('centroid')

        hold on
%         if maptype == 1
%         scatter(temp_data(:,3)-temp_data(:,1),temp_data(:,4)-temp_data(:,2),5,ori_pref,'')
%         end
                axis equal
        set(gca,'YDir','reverse')
    end
   

end
