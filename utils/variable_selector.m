function [target_data] = variable_selector(targets,str,selection_vector)

% separate the targets
targets = strsplit(targets,',');

% get the number of targets
number_targets = numel(targets);

% allocate memory for the target data
target_data = zeros(sum(selection_vector),number_targets);
% for all the targets
for target = 1:number_targets

    if contains(targets{target},'ang_')

        % parse the argument
        field_name = targets{target}(1:end-3);
        column_name = targets{target}(end-1:end);
        % load the field
        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'Cx'
                target_data(:,target) = abs(temp_var(:,3)-temp_var(:,1));
            case 'Cy'
                target_data(:,target) = abs(temp_var(:,4)-temp_var(:,2));
            case 'Vt'
                target_data(:,target) = temp_var(:,8);
            case 'Al'
                target_data(:,target) = temp_var(:,5);
            case 'Sx'
                target_data(:,target) = temp_var(:,1);
            case 'Sy'
                target_data(:,target) = temp_var(:,2);
            case 'Rx'
                target_data(:,target) = temp_var(:,3);
            case 'Ry'
                target_data(:,target) = temp_var(:,4);
            case 'Dx'
                target_data(:,target) = temp_var(:,3)-temp_var(:,1);
            case 'Dy'
                target_data(:,target) = temp_var(:,4)-temp_var(:,2);
        end
    elseif contains(targets{target},'nw_')
        % parse the argument
        field_name = strcat('ang',targets{target}(3:end-3),'_nonweighted');
        column_name = targets{target}(end-1:end);
        % load the field
        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'Cx'
                target_data(:,target) = abs(temp_var(:,3)-temp_var(:,1));
            case 'Cy'
                target_data(:,target) = abs(temp_var(:,4)-temp_var(:,2));
            case 'Vt'
                target_data(:,target) = temp_var(:,8);
            case 'Al'
                target_data(:,target) = temp_var(:,5);
            case 'Sx'
                target_data(:,target) = temp_var(:,1);
            case 'Sy'
                target_data(:,target) = temp_var(:,2);
            case 'Rx'
                target_data(:,target) = temp_var(:,3);
            case 'Ry'
                target_data(:,target) = temp_var(:,4);
            case 'Dx'
                target_data(:,target) = temp_var(:,3)-temp_var(:,1);
            case 'Dy'
                target_data(:,target) = temp_var(:,4)-temp_var(:,2);
        end
    elseif contains(targets{target},{'frac_vert','frac_horz'})
        % parse the argument
        field_name = targets{target}(1:9);
        column_name = targets{target}(11:end);

        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'exL23'
                target_data(:,target) = sum(temp_var(:,3:5),2);
            case 'exL4'
                target_data(:,target) = sum(temp_var(:,6:7),2);
            case 'exL5'
                target_data(:,target) = sum(temp_var(:,8:10),2);
            case 'inL23'
                target_data(:,target) = sum(temp_var(:,19:21),2);
            case 'inL4'
                target_data(:,target) = sum(temp_var(:,22:23),2);
            case 'inL5'
                target_data(:,target) = sum(temp_var(:,24:26),2);
        end
    elseif contains(targets{target},'custom')
        temp_var = vertcat(str(selection_vector).ang_inL23);
        x = abs(temp_var(:,3)-temp_var(:,1));
        y = abs(temp_var(:,4)-temp_var(:,2));
        target_data = sqrt(x.^2 + y.^2);
    elseif contains(targets{target},'max_')
        % parse the argument
        field_name = targets{target}(1:6);
        coord_name = targets{target}(8);
        layer_name = targets{target}(10:end);

        % get the data
        temp_var = reshape(vertcat(str(selection_vector).(field_name)),[],3,3);
        switch coord_name
            case 'v'
                coord = 1;
            case 'x'
                coord = 2;
            case 'y'
                coord = 3;
        end
        switch layer_name
            case 'L23'
                layer = 1;
            case 'L4'
                layer = 2;
            case 'L5'
                layer = 3;
        end

        % store the data
        target_data(:,target) = squeeze(temp_var(:,coord,layer));
    elseif contains(targets{target},'span')
        % parse the argument
        parts = strsplit(targets{target},'_');
        layer_name = parts{2};
        switch layer_name
            case 'L23'
                layer_coord = 1;
            case 'L4'
                layer_coord = 2;
            case 'L5'
                layer_coord = 3;
        end
        % if it's inhibition, offset the coordinate
        if strcmp(parts{3},'in')
            layer_coord = layer_coord + 3;
        end 
        % load the data
        temp_data = vertcat(str(selection_vector).span);
        % output the correct column
        target_data(:,target) = temp_data(:,layer_coord);

    else
        target_data(:,target) = vertcat(str(selection_vector).(targets{target}));
    end
end