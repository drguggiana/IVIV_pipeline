% Quantify the overlap between input map and neuronal morphology
%% Clean up
clearvars
close all
%% Load the neuronal trees

% % define the path for the trees
% tree_path = 'R:\Share\Simon\Drago_Volker_Simon\InVivivo_InVitro full data structure\190612_1400_dataStruct.mat';
% % define the save path
% save_path = 'R:\Share\Simon\Drago_Volker_Simon\morphoMaps\morphoMaps.mat';
Paths
% load the structure
str = load(structure_file_path);
str = str.str;
%% Grid definition

% % define the grid size
% grid_size = 16;
% % define the grid spacing in microns
% grid_spacing = 69;
% % generate the grid centered on zero
% x_vector = linspace(0, grid_spacing*grid_size, grid_size);
% x_vector = x_vector - max(x_vector/2);
% y_vector = x_vector;
% 
% % generate the points
% [X,Y] = meshgrid(x_vector,y_vector);
% paired_points = cat(2,X(:),Y(:));
% 
% % load the xsg info
% xsg_info = load('R:\Share\Simon\Drago_Volker_Simon\xsg_info.mat');
% xsg_info = xsg_info.xsg_info;
%% Generate image of a neuron
close all

% define the target interpolations density [pt/px]
target_density = 10;

% define the neuron indexes to work with
% neuron_idx = 1:10;
neuron_idx = 1:length(str);
neuron_num = length(neuron_idx);
% allocate a structure to save the data
interp_neurons = zeros(neuron_num,1);

% define the histogram bins
% bins = [0 45 90 135 180 225 270 315 360 405]-45/2;
% bins = [0 45 90 135 180 225]-45/2;
bins = linspace(-100,100,10);

% allocate memory for the histograms
histogram_matrix = cell(neuron_num,1);
% allocate memory for the angles
angle_cell = cell(neuron_num,1);
% generate and image for the neuron
% for all neurons in the vector
for neuron = 1:130%neuron_num    
    %%
    % get the index of the neuron
    idx = neuron_idx(neuron);
    
    % get the tree info
    curr_neuron = str(idx).morphtraces;
    % if there is no morpho, skip the cell
    if isempty(curr_neuron) || ~iscell(curr_neuron) || all(isnan(str(idx).morph)) || neuron==147
        continue
    end
    
    ori_pref = str(neuron).ORIpref;
    % if it's NaN, skip
    if isnan(ori_pref)
        continue
    end
    % allocate memory for the neurites
    neurite_cell = cell(3,1);
    % for all neurite types
    for neurite = 1:2
        %%
        % get the neurite
        curr_neurite = curr_neuron{neurite};
        
        % flip it on the y axis
        curr_neurite.Y = -curr_neurite.Y;
        
        % if the slice is inverted, flip the cell
        if str(neuron).sliceOri == 1
            curr_neurite.X = -curr_neurite.X;
            curr_neurite.Z = -curr_neurite.Z;
        end
        
        % store the center of mass of the cell
        if neurite == 3
            cell_center = [mean(curr_neurite.X), mean(curr_neurite.Y)];
        end
        % plot a single neurite
%         figure
%         plot3(curr_neurite.X, curr_neurite.Y, curr_neurite.Z, 'o')
        % interpolate the neurite segments to a common density
%         figure
%         imagesc(curr_neurite.dA)
        % get the node indexes
        [node_1, node_2] = find(curr_neurite.dA);
        node_number = length(node_1);
        % start a counter for the each of the nodes
        node_counter = zeros(node_number+1,1);
        distance = zeros(node_number,1);
        % allocate memory for the interpolated segments
        segment_cell = cell(node_number, 1);
%         figure
        % for all the nodes
        for nodes = 1:node_number
            % count the start
            node_counter(node_1(nodes)) = node_counter(node_1(nodes)) + 1;
           
            % get the node coordinates for start and end
            node_start = [curr_neurite.X(node_1(nodes)), ...
                curr_neurite.Y(node_1(nodes)), ...
                curr_neurite.Z(node_1(nodes))];
            node_end = [curr_neurite.X(node_2(nodes)), ...
                curr_neurite.Y(node_2(nodes)), ...
                curr_neurite.Z(node_2(nodes))];
            % get the distance between the start and end
            node_distance = norm(node_end - node_start);
            distance(nodes) = node_distance;
            % if the counter value is > 2 or the distance is too long, skip
            if node_counter(node_1(nodes)) > 2 || node_distance > 20
                continue
            end
            
            % calculate the number of points for the segment
            interp_points = round(node_distance*target_density);
            % assemble the segment
            segment_cell{nodes} = [linspace(node_start(1), node_end(1), interp_points);...
                linspace(node_start(2), node_end(2), interp_points); ...
                linspace(node_start(3), node_end(3), interp_points)];
            
            
%             plot3(segment_cell{nodes}(1,:),segment_cell{nodes}(2,:),segment_cell{nodes}(3,:),'o')
%             hold('on')
            
        end
        
        % save the segments
%         neurite_cell{neurite} = cat(2,segment_cell{:});
        neurite_cell{neurite} = [curr_neurite.X,curr_neurite.Y,curr_neurite.Z]';
        
%         close all
%         
%         scatter3(neurite_cell{neurite}(1,:),neurite_cell{neurite}(2,:),neurite_cell{neurite}(3,:))
%         figure
%         histogram(distance)
    end
    % combine the neuron coordinates
    single_neuron = cat(2,neurite_cell{:});
    
    % trim to a quarter of the neuron
%     single_neuron = single_neuron(:,single_neuron(1,:)>0);
%     single_neuron = single_neuron(:,single_neuron(3,:)>0);
    
%     close all
%     figure
%     scatter3(single_neuron(1,:),single_neuron(2,:),single_neuron(3,:))
%     title(num2str(neuron))
%     axis equal
%     scatter(single_neuron(1,:),single_neuron(2,:),5)
%     title(num2str(neuron))
%     axis equal
%     hold on
%     ellipse = fit_ellipse(single_neuron(1,:),single_neuron(3,:));
%     
%     % save the value
%     interp_neurons(neuron) = rad2deg(ellipse.phi)+90;
    
    %% calculate radial density
    
    % get the angle of every point
    [angles,rho] = cart2pol(single_neuron(1,:),single_neuron(3,:));
    
    % define the rho threshold
    rho_threshold = 20;
    angles = rad2deg(angles(rho<rho_threshold));
    
%     figure
%     scatter(single_neuron(1,rho<rho_threshold),single_neuron(3,rho<rho_threshold),5)
%     title(num2str(neuron))
%     axis equal
%     hold on

    figure
    scatter3(single_neuron(1,rho<rho_threshold),single_neuron(2,rho<rho_threshold),single_neuron(3,rho<rho_threshold))
    title(num2str(neuron))
    axis equal
    
    % exclude angles beyond a radius
    
    % rotate by the orientation preference of the cell
    ori_pref = str(neuron).ORIpref;
    ori_pref = ori_pref + 90;
    ori_pref(ori_pref>179.999999) = ori_pref(ori_pref>179.999999) - 180;

    
    % otherwise rotate and wrap
%     angles = angles - ori_pref;
%     angles(angles>180) = angles(angles>180) - 180;
%     angles(angles<0) = angles(angles<0) + 180;

%     angles(angles>90) = angles(angles>90) - 180;
%     angles(angles<-90) = angles(angles<-90) + 180;


    % store the angles
    angle_cell{neuron} = angles;
    
    % get the histogram and store
    N = histcounts(angles, bins);
    
%     % combine the first and last bins
%     N = [N(2:end-1), N(1)+N(end)];
    
    % store
    histogram_matrix{neuron} = N;
    
end
% concatenate the histogram matrix
histogram_matrix = vertcat(histogram_matrix{:});
% concatenate the angles
angle_cell = horzcat(angle_cell{:});
autoArrangeFigures
%% Plot the histogram

close all

% figure
% imagesc(histogram_matrix)

figure
errorbar(bins(1:end-1)+5, mean(histogram_matrix,1), std(histogram_matrix,0,1)./sqrt(size(histogram_matrix,1)))

figure
polarhistogram(deg2rad(angle_cell+180))

autoArrangeFigures
%% Compare to ori pref

close all
cell_selection = [str.OSIpref]>0.25 & abs(interp_neurons')>0;

figure
scatter([str(cell_selection).ORIpref], interp_neurons(cell_selection))

[rho,pval] = circ_corrcc(deg2rad([str(cell_selection).ORIpref]), deg2rad(interp_neurons(cell_selection)))
%% Save the output as a file
% save(save_path,'interp_neurons')