%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = old_structure_file_path;
str = load(main_path);
str = str.str;
% get the number of cells
cell_num = size(str,1);
%% Get rid of the 9 cells with NaN inhibition maps

% allocate a vector to store the cells to keep
keep_cells = ones(cell_num,1);
% for all the cells
for cells = 1:cell_num
    % check if either the inhibiton or the excitation maps are NaN, if so,
    % record it
    if isnan(sum(str(cells).excMap(:))) || isnan(sum(str(cells).inhMap(:)))
        keep_cells(cells) = 0;
    end
end

% exclude the cells
str = str(keep_cells==1);

% update the cell_num
 cell_num = size(str,1);
%% Morpho density corr
% calculate the correlation between morpho density and maps

% % allocate memory for the correlations (exc/inh,apical/basal)
% correlation_values = zeros(cell_num,2,2);

% for all the cells
for cells = 1:cell_num

    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = str(cells).excMap*-1;
                term1 = 'exc';
            case 2
                map = str(cells).inhMap;
                term1 = 'inh';
        end
        % get rid of layer 1
        map = map(3:end,:);
        % for both types of morpho map
        for morph = 1:2
            switch morph
                case 1
                    morpho = str(cells).morphoMap_apical_aligned;
                    term2 = 'apical';
                case 2
                    morpho = str(cells).morphoMap_basal_aligned;
                    term2 = 'basal';
            end
            % if morpho is empty, skip and load a NaN
            if isempty(morpho)
                str(cells).(strjoin({term1,term2},'_')) = NaN;
                continue
            end
            morpho = morpho(3:end,:);
            str(cells).(strjoin({'corr',term1,term2},'_')) = corr(map(:), morpho(:));
        end
        
    end
end
%% Add the noise correlations

% get the list of folders
folder_list_raw = dir(od_svd_path);
%get the number of folders
exp_num = length(folder_list_raw);
% allocate memory for the full paths
folder_list = cell(exp_num,1);
% for all the folders
for folders = 1:exp_num
    folder_list{folders} = fullfile(folder_list_raw(folders).folder,...
        folder_list_raw(folders).name);
end
%allocate memory to store the ROI and neuropil data
svd_OD = cell(exp_num,1);
% also for the noise correlations
noise_OD = cell(exp_num,1);

% allocate memory to also save the number of rois per experiment along with
% the experiment name
cell_number = cell(exp_num,2);

%for all the folders
for experiment = 1:exp_num
    % load the cell containing the hosvd decompositions
    svd_OD{experiment} = load(folder_list{experiment});
    noise_OD{experiment} = svd_OD{experiment}.noise_matrix;
    % collapse the second cell dimension
    noise_OD{experiment} = squeeze(mean(noise_OD{experiment},2));
    svd_OD{experiment} = svd_OD{experiment}.cell_cell;
    
    [~,cell_number{experiment,1}] = fileparts(folder_list{experiment});
    cell_number{experiment,1} = cell_number{experiment,1}(4:8);
    cell_number{experiment,2} = size(svd_OD{experiment},1);
end

% concatenate the data
svd_OD = cat(1,svd_OD{:});
noise_OD = cat(1,noise_OD{:});
%% Match the iviv cells with the all invivo cells

% load the excel spreadsheet with the matching
[~,~,matching_raw] = xlsread(matching_path);

% get the animal names and their corresponding OD names
matching_names = matching_raw(2:end,[4,9,11]);
% % turn the experiment strings into numbers
matching_names(:,2) = cellfun(@eval, matching_names(:,2),'UniformOutput',0);
% for all the rows (cause NaNs in the cellID field -_-
for rows = 1:size(matching_names,1)
    if ~ischar(matching_names{rows,3})
        matching_names{rows,3} = num2str(matching_names{rows,3});
    end
    matching_names{rows,3} = eval(matching_names{rows,3});
end

% allocate memory to store the pial depth
iviv_vector = cell(exp_num,1);

% for all the experiments
for experiment = 1:exp_num
    % get the current experiment name
    [~,current] = fileparts(folder_list{experiment});
    % find the matching animal name from the iviv names
    iviv_idx = contains(matching_names(:,1),current(4:8));
    % assemble a vector with the cell number
    iviv_temp = zeros(cell_number{contains(cell_number(:,1),current(4:8)),2},2);
    % turn the iviv indexes into ones
    iviv_temp(matching_names{iviv_idx,2},1) = 1;
    % and the cellID fields into their corresponding number (only if not
    % NaN)
    if ~isnan(matching_names{iviv_idx,3})
        iviv_temp(matching_names{iviv_idx,2},2) = matching_names{iviv_idx,3};
    end
    % store the resulting vector in the main cell
    iviv_vector{experiment} = iviv_temp;
end

% concatenate the cell
iviv_vector = cat(1,iviv_vector{:});
% % get the corrected cells only
% iviv_vector = iviv_vector(correct_cells,:);
% iviv_vector = iviv_vector(correct_cells_SFTF,:);
%% Save the noise correlations

% get the noise correlations only for the matched cells
% noise_matched = noise_OD(correct_cells,:,:);
% noise_matched = noise_matched(correct_cells_SFTF,:,:);
noise_matched = noise_OD(iviv_vector(:,1)==1,:,:);

% get the cell IDs
cell_id = cat(1,str.cellID);
% get ids from the iviv vector
iviv_id = iviv_vector(iviv_vector(:,1)==1,2);

% for all the matched cells
for cells = 1:size(noise_matched,1)
    % get the boolean for selection of the cell
    id_bool = cell_id==iviv_id(cells);
    % if the cell is not here, print the id and skip
    if sum(id_bool)==0
        fprintf(strjoin({'Cell absent:',num2str(iviv_id(cells)),'\r\n'},'_'))
        continue
    end
    % take the correlation for the preferred direction
    [~,~,bin] = histcounts(str(id_bool).DIRpref,-22.5:45:382.5);
    % rectify the last bin
    if bin == 9
        bin = 1;
    elseif bin == 0
        str(id_bool).noise = NaN;
        continue
    end
    % check the preference of the cell and take the corresponding average
    if str(id_bool).contra == 1
        str(id_bool).noise = squeeze(noise_matched(cells,bin,1));
    elseif str(id_bool).ipsi == 1
        str(id_bool).noise = squeeze(noise_matched(cells,bin,2));
    else
        if str(id_bool).ODIpref > 0
            str(id_bool).noise = squeeze(noise_matched(cells,bin,1));
        else
            str(id_bool).noise = squeeze(noise_matched(cells,bin,2));
        end
    end
end
%% Turn all the empty fields into NaN

% get the fields
field_list = fields(str);
% get the number of fields
num_fields = length(field_list);
% for all the fields
for f = 1:num_fields
    % for all the cells
    for cells = 1:cell_num
        if isempty(str(cells).(field_list{f}))
            str(cells).(field_list{f}) = NaN;
        end
    end
end
%% Save the structure (keep moving this down as features are added)
save(structure_file_path,'str')