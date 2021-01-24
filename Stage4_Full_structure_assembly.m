%% Collapse the in vivo, in vitro and morphology structures in one file

%% Clean up 

clearvars
close all

% get the paths
Paths
%% Load the component structures

% load the in vitro structure
invitro_struct = load(find_newer_file(stage1_invitro_path));
invitro_struct = invitro_struct.invitro_struct;

% load the in vivo structure
invivo_struct = load(find_newer_file(stage2_invivo_path));
invivo_struct = invivo_struct.invivo_struct;

% load the mephys struct
mephys_struct = load(find_newer_file(stage3_mephys_path));
mephys_struct = mephys_struct.mephys;
%% Match the cells

% get the cell names
invitro_names = {invitro_struct.cellName}';
invivo_names = {invivo_struct.cellnames}';
invivo_names = vertcat(invivo_names{:});
mephys_names = {mephys_struct.cell_name}';

% get the total number of cells (i.e. unique)
overall_names = unique(vertcat(invitro_names,invivo_names,mephys_names));

% get the total number of cells
total_cells = length(overall_names);
%% Assemble the overall structure

% allocate the structure
str_allcells = struct([]);

% get a list of the fields for each structure
invitro_fields = fields(invitro_struct);
invitro_fields = invitro_fields(~contains(invitro_fields,{'cellName'}));
invivo_fields = fields(invivo_struct);
invivo_fields = invivo_fields(~contains(invivo_fields,{'cellnames'}));
mephys_fields = fields(mephys_struct);
mephys_fields = mephys_fields(~contains(mephys_fields,{'cell_name'}));

% for all the cells
for cells = 1:total_cells
    % get the cell name
    current_cell = overall_names{cells};
    
    % save the cell name
    str_allcells(cells).cellName = current_cell;
    % initialize the cellID field
    str_allcells(cells).cellID = [];
    
    % if the cells is in vitro, copy fields
    if any(contains(invitro_names,current_cell))
        % get the cell's invitro idx
        idx = find(contains(invitro_names,current_cell));
        % for all the in vitro fields
        for field = 1:length(invitro_fields)
            % copy the info
            str_allcells(cells).(invitro_fields{field}) = ...
                invitro_struct(idx).(invitro_fields{field});
        end
    end
    
    % if the cells is in vivo, copy fields
    if any(contains(invivo_names,current_cell))
        % get the cell's invivo idx
        idx_invivo = find(contains(invivo_names,current_cell));
        % for all the in vivo fields
        for field = 1:length(invivo_fields)
            % if the field is the id, overwrite the cellID field from
            % invitro (should be the same number)
            if strcmp(invivo_fields{field},'cellids')
                
                str_allcells(cells).cellID = ...
                    invivo_struct(idx_invivo).(invivo_fields{field});
            else
                % copy the info
                str_allcells(cells).(invivo_fields{field}) = ...
                    invivo_struct(idx_invivo).(invivo_fields{field});
            end
        end
    else
        % set as empty to identify it as non matching below
        idx_invivo = [];
    end
    
    % if the cells is in vitro, copy fields
    if any(contains(mephys_names,current_cell))
        % get the cell's mephys idx
        idx = find(contains(mephys_names,current_cell));
        % for all the in vitro fields
        for field = 1:length(mephys_fields)
            % if the field is the id, overwrite the cellID field from
            % invitro (should be the same number)
            if strcmp(mephys_fields{field},'id')
                str_allcells(cells).cellID = ...
                    mephys_struct(idx).(mephys_fields{field});
            else
                % copy the info
                str_allcells(cells).(mephys_fields{field}) = ...
                    mephys_struct(idx).(mephys_fields{field});
            end
        end
    end
end
%% Save the structure

% define the save path
save_path = fullfile(stage4_full_structure_path,'all_cells.mat');

% save the file
save(save_path,'str_allcells');