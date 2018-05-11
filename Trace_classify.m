%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load paths

%pick the folders to use
folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\INPUT MAPS_final\');


%define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';
%get the number of folders
folder_num = length(folder_list);
%allocate memory to store the subfolders that qualify
folder_cell = cell(folder_num,1);
%get the setup of origin (for later)
ori_setup = strsplit(folder_list{1},'\');
ori_setup = ori_setup{end-1};

%for all the folders
for folders = 1:folder_num

    %get the subfolders in the first level of this folder
    subfolder_list = dir(folder_list{folders});

    %also, get rid of the non-folders
    subfolder_list = subfolder_list(vertcat(subfolder_list(:).isdir)==1);
    %leave only the folders that have SW in them
    subfolder_list = subfolder_list(contains({subfolder_list(:).name},{'SW','FM'}));
    %finally, build the full paths for the folders, including the "images"
    %folder and the "map01"
    path_list = cell(size(subfolder_list,1),1);
    %for all the subfolders
    for subs = 1:size(subfolder_list,1)

        %load all the map paths in a given cell
        
        %get the list of paths
        xsg_paths = dir(strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\*map*'));
        %and create a cell to store the full paths
        xsg_cell = cell(length(xsg_paths),1);
        %also a vector to keep track of the empty folders
        elim_vec = ones(size(xsg_cell,1),1);
        %for each path
        for paths = 1:length(xsg_paths)
            %load the files in the directory
            xsg_files = dir(strcat(xsg_paths(paths).folder,'\',xsg_paths(paths).name,'\*.xsg'));
            %if there is no file, skip the entry and update the vector
            %accordingly
            if isempty(xsg_files)
                elim_vec(paths) = 0;
                continue
            end
            %append the name of the file to the path
            xsg_cell{paths} = strcat(xsg_files(1).folder,'\',xsg_files(1).name);
        end
        %store the cell with paths in the larger storage cell
        path_list{subs} = xsg_cell(elim_vec == 1);
    end
%     %update the folder cell only with the cells that had an image
%     folder_cell{folders} = path_list(elim_vec==1,:);
    %update the folder cell
    folder_cell{folders} = vertcat(path_list{:});
end

%concatenate the entire list of paths
folder_all = vertcat(folder_cell{:});
%get the number of maps to do
map_num = size(folder_all,1);
%% Get an array with the list of cell names for each map

%extract the cell names from the folder all array
%allocate memory to save the cell names
cell_names = cell(map_num,1);
%for all the folders
for maps = 1:map_num
    temp_name = strsplit(folder_all{maps},'\');
    cell_names{maps} = strcat(temp_name{end-2}(1:2),temp_name{end-3}(1:6),temp_name{end-2}(3:6));
end
%% Load, trim and reorder the data for the PCA 

%load a random 10% of the data and run PCA on it
%get the random indexes of the maps to load
% rand_ind = randperm(map_num,round(map_num/10));
% rand_ind = randperm(map_num,68);
%allocate memory to store the maps
temp_xsg = load(folder_all{1},'-mat');
% %get the size of the maps
% map_size = length(temp_xsg.data.ephys.trace_1);
%also the number of positions
num_positions = temp_xsg.header.mapper.mapper.positionNumber;
%define the interval to grab. Each recording is a second and the desired
%range is 150 ms after the direct window (which occurs 7 ms after the 100ms
%post stimulus onset)
trace_range = 1001:3000;
%define the time interval for the background
trace_background = 1:1000;
%define the time interval for the direct window
trace_direct = 1001:1070;
trace_pseudodirect = 1001:1036;
%allocate a matrix for them
% map_matrix = zeros(length(trace_range),length(rand_ind));
map_matrix = zeros(length(trace_range)*num_positions,map_num);
%allocate memory for the standard deviation of the background for each
%trace
background_std = zeros(num_positions,map_num);
%allocate memory for the direct and pseudodirect window std
direct_std = zeros(num_positions,map_num);
pseudodirect_std = zeros(num_positions,map_num);
%allocate memory for the soma centers
soma_center = zeros(map_num,2);

%set up a vector for the skipped files
elim_vec = ones(size(map_matrix,2),1);

%create a waitbar
w_bar = waitbar(0,'Loading maps');
%load the random maps into the matrix
% for maps = 1:length(rand_ind)
for maps = 1:map_num
    
    %show progress
    waitbar(maps/map_num,w_bar)
    %load the xsg file
%     temp_xsg = load(folder_all{rand_ind(maps)},'-mat');
    temp_xsg = load(folder_all{maps},'-mat');
    %check the number of grid points in the cell
    temp_map = temp_xsg.header.mapper.mapper.mapPatternArray;
    %if the trace has less than 16x16
    if size(temp_map,1) < 16 || size(temp_map,2) < 16 ||...
            isempty(temp_xsg.data.ephys)
        %mark the position and skip the trace
        elim_vec(maps) = 0;
        continue
    end
    %load the time trace
    temp_trace = temp_xsg.data.ephys.trace_1;
    %reshape to trim the relevant part of the data only
    temp_trace = reshape(temp_trace,[],num_positions);
   
    %get the background activity
    background_act = mean(temp_trace(trace_background,:),1);
    %get the std deviation of the background
    background_std(:,maps) = std(temp_trace(trace_background,:),0,1);
    %calculate the direct window std
    direct_std(:,maps) = std(temp_trace(trace_direct,:),0,1);
    pseudodirect_std(:,maps) = std(temp_trace(trace_pseudodirect,:),0,1);
    %trim the trace
    temp_trace = temp_trace(trace_range,:);
    %subtract background activity
    temp_trace = temp_trace - background_act;
    %also load the map order to have all the maps ordered in their location
    %instead of in time of stimulation
    %linearize the map and apply to the data
    temp_trace = temp_trace(:,temp_map(:));
    %reorder the background and direct maps too 
    background_std(:,maps) = background_std(temp_map(:),maps);
    direct_std(:,maps) = direct_std(temp_map(:),maps);
    pseudodirect_std(:,maps) = pseudodirect_std(temp_map(:),maps);
    %reshape again and store
    map_matrix(:,maps) = reshape(temp_trace,length(trace_range)*num_positions,[]);
    %if not empty,save the soma center(assume that at least one of the xsg
    %files for the cells will have it)
    if ~isempty(temp_xsg.header.mapper.mapper.soma1Coordinates)
        soma_center(maps,:) = temp_xsg.header.mapper.mapper.soma1Coordinates;
    else
        soma_center(maps,:) = NaN;
    end
end

close(w_bar)

%get rid of the empty spaces with the skipped traces
map_matrix = map_matrix(:,elim_vec==1);
%also within the std matrix
background_std = background_std(:,elim_vec==1);
direct_std = direct_std(:,elim_vec==1);
pseudodirect_std = pseudodirect_std(:,elim_vec==1);
% map_matrix2 = map_matrix;
%and also modify the folder vector (so I can refer back to the particular
%maps)
folder_all = folder_all(elim_vec==1);
%update the cell name array
cell_names = cell_names(elim_vec==1);
%also the soma centers
soma_center = soma_center(elim_vec==1,:);

%get the unique cell names
[uni_cells,ia,cell_id] = unique(cell_names);
%replace the empty centers by the next value (assuming the first file was
%messed up only)
% empty_idx = find(isnan(soma_center(:,1)));
% soma_center(empty_idx,:) = soma_center(empty_idx+1,:);
%get the unique soma centers
soma_unique = soma_center(ia,:);
%% Split the maps into excitatory and inhibitory

%define the map path
excinh_path = 'R:\Share\Simon\Drago_Volker_Simon\';

%load the file
load(fullfile(excinh_path,strcat('ExcInh_map_',ori_setup,'.mat')))

%get only the binary values for the responses
map_polarity = vertcat(ExcInh_map{:,1});

% %calculate the average value for each map, since area under the curve
% %should report whether it is exc or inh
% ave_value = mean(map_matrix,1);
% 
% figure
% histogram(ave_value)
% %calculate the polarity of each map to propagate later (0 is exc)
% map_polarity = ave_value>0;
%% Reshape the matrix, so that all traces are in a single column

map_matrix2 = reshape(map_matrix,length(trace_range),[]);
%linearize the std matrices also
background_std = background_std(:);

%and create a map from traces to folders (including polarity)
%allocate memory for the map
trace2folder = zeros(num_positions,size(map_matrix,2),6);
%for all the folders
for folders = 1:size(folder_all,1)
    trace2folder(:,folders,1) = 1:num_positions;
    trace2folder(:,folders,2) = folders;
    trace2folder(:,folders,3) = map_polarity(folders);
    trace2folder(:,folders,4) = cell_id(folders);
    
end
%also reshape this map to be able to refer to the original map
trace2folder = reshape(trace2folder,[],6);
%% Filter out traces that are too flat (using 3 std criterion)
close all
%calculate the std of each trace
std_matrix = std(map_matrix2,0,1);

%plot the distribution of the std
figure
histogram(std_matrix)

%get the direct and direct mix traces
direct_std = direct_std(:);
direct_map = direct_std>3.*background_std;

pseudodirect_std = pseudodirect_std(:);
pseudodirect_map = pseudodirect_std>3.*background_std;

%Calculate the comparison between negative and positive deflections in the
%inh traces (1 is positive deflection is higher)
% deflection_map = (max(map_matrix2,[],1)>abs(min(map_matrix2,[],1)))';
deflection_map = zeros(size(map_matrix2,2),1);
%for all the traces
for traces = 1:size(map_matrix2,2)
    deflection_map(traces) = sum(map_matrix2(map_matrix2(:,traces)>0,traces))...
        >-sum(map_matrix2(map_matrix2(:,traces)<0,traces));
end
%and store the direct list also
trace2folder(:,6) = direct_map;
% %exclude traces with less than 3 stds over background
% map_matrix2 = map_matrix2(:,std_matrix'>3.*background_std);
%also save a map of the selected traces
prctile_map = std_matrix'>3.*background_std;
% %generate a map of the origin of the traces left after the percentile cut
% postperc_map = trace2folder(prctile_map,:);
%mark the places on the trace2folder matrix
trace2folder(prctile_map==0,5) = NaN;
%% Fill in the classification vector depending on the window criteria


%for both polarities
for polarity = 1:2
    
    %get the idx from this polarity
    polarity_idx = trace2folder(:,3)==polarity-1;
    
    %switch depending on inh vs exc traces (in case we decide to do
    %something different for the two)
    switch polarity
        case 1 %excitatory maps
            %mark the pure synaptic traces on the 5th column
            trace2folder(prctile_map==1&polarity_idx&direct_map==0,5) = 1;
            
            %now mark the pseudodirect ones in the same column (leaving true synaptic
            %as 1 and pseudo as 2)
            trace2folder(prctile_map==1&polarity_idx&direct_map==1&pseudodirect_map==0,5) = 2;
            
        case 2 %inhibitory maps
            
%             %mark the pure synaptic traces on the 5th column
%             trace2folder(prctile_map==1&polarity_idx&direct_map==0,5) = 1;
% 
%             %now mark the pseudodirect ones in the same column (leaving true synaptic
%             %as 1 and pseudo as 2)
%             trace2folder(prctile_map==1&polarity_idx&direct_map==1&pseudodirect_map==0,5) = 2;
            %mark the traces that have a higher positive deflection, syn
            %window ok
            trace2folder(prctile_map==1&polarity_idx&deflection_map==1&direct_map==0,5) = 1;
            
            %mark the traces where the negative deflection is higher, syn
            %window ok
            trace2folder(prctile_map==1&polarity_idx&deflection_map==0&direct_map==0,5) = 2;
            
            %mark the traces that have a higher positive deflection, syn
            %window out
            trace2folder(prctile_map==1&polarity_idx&deflection_map==1&direct_map==1,5) = 3;
            
    end

end
%% Save the results

% %define the save path
% save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_rawClass',ori_setup,'.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(out_path,save_name),'folder_all','prctile_map',...
    'trace_background','trace_range','uni_cells',...
    'trace2folder','background_std','trace_direct','trace_pseudodirect','soma_unique')
% error('Done')