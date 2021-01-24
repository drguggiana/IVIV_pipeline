%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;

% load the invivo structure
invivo = load(invivo_path);
invivo = invivo.L23_PC;

% load the correction factor for the pial depth
pia_correction = load(pia_correction_path);
pia_correction = pia_correction.scale_di;
%% Parse the invivo structure for what I need

% get the number of experiments
experiment_number = size(invivo,2);

% define the target list of fields
field_list = {'prefOri','prefDir','gOSI','gDSI','ODI'};
invivo_str = struct([]);
cell_counter = 1;
% for all the experiments
for experiment = 1:experiment_number

    % get the current structure
    current_str = invivo(experiment).OD;
    % get the pia
    current_pia = invivo(experiment).pial_depth;
    % correct based on the correction factor
    current_pia = current_pia - pia_correction(experiment);
    % get the number of cells
    cell_num = size(current_str.contra,2);
    % start a counter

    % for all the cells
    for cells = 1:cell_num
        % if it's a non-responder, skip
        if current_str.oresp(cells) == 0
            continue
        end
        % determine the preferred parameters
        if current_str.contra(cells) == 1
            selector = 1;
        elseif current_str.ipsi(cells) == 1
            selector = 2;
        else
            if current_str.ODI(cells) > 0
                selector = 1;
            else
                selector = 2;
            end
        end
        
        % for all the fields
        for field = 1:size(field_list,2)
            % get the field
            current_field = current_str.(field_list{field});
            % if there are 2 dimensions, use the pref one
            if size(current_field,2)==2
                current_field = current_field(cells,selector);
            else
                current_field = current_field(cells);
            end
            % save the info
            invivo_str(cell_counter).(field_list{field}) = current_field;
        end
        % also add the pia
        invivo_str(cell_counter).pialD = current_pia(cells);
        % update the counter
        cell_counter = cell_counter + 1;
    end
end
%% Pia plot
close all

% correct cell number
num_cells = size(invivo_str,1);

% get the OSI
OSI = vertcat(invivo_str.gOSI);

% get the orientation preference
ORI_pref = vertcat(invivo_str(OSI>0.25).prefOri);

% turn the values into bar orientation
% ORI_pref(ORI_pref>180) = ORI_pref(ORI_pref>180) - 180;
ORI_pref = ORI_pref + 90;
ORI_pref(ORI_pref>179.999999) = ORI_pref(ORI_pref>179.999999) - 180;

% get the pia and also filter
piald = vertcat(invivo_str(OSI>0.25).pialD);

% separate the orientation in bins based on pia
% [N,edges,bin] = histcounts(piald,[100 200 300 500]);
[N,edges,bin] = histcounts(piald,[100 200 250 500]);
% [N,edges,bin] = histcounts(piald,3);

% define edges for the angle bins
% angle_edges = [0 45 90 135 180 225]-22.5;
% angle_edges = [0 32.5 82.5 122.5 172.5 225]-22.5;
angle_edges = [0 35 80 125 170 225]-22.5;
% allocate memory for the bin results
depth_bins = zeros(size(N,2),size(angle_edges,2)-2);
% for all the depth bins
for bins = 1:size(N,2)
    % bin the angles
    temp_bins = histcounts(ORI_pref(bins==bin),angle_edges,'Normalization','probability');
    % combine the edge bins
    depth_bins(bins,:) = [temp_bins(2:4),sum(temp_bins([1,5]))];
end

% plot the results
bar(depth_bins')
set(gca,'XTickLabels',string(angle_edges(2:5)+22.5))
legend({'High','Mid','Low'})
xlabel('Preferred Orientation')
ylabel('Probability')