%% Clean up
clearvars
close all
addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load paths

%pick the folders to use
folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\INPUT MAPS_final\Setup2');


% %define the path to save the output files
% out_path = 'R:\Share\Simon\Drago_Volker_Simon\layer_GUI_out';
%get the number of folders
folder_num = length(folder_list);
%allocate memory to store the subfolders that qualify
folder_cell = cell(folder_num,1);
%for all the folders
for folders = 1:folder_num

    %get the subfolders in the first level of this folder
    subfolder_list = dir(folder_list{folders});

    %also, get rid of the non-folders
    subfolder_list = subfolder_list(vertcat(subfolder_list(:).isdir)==1);
    %leave only the folders that have SW in them
    subfolder_list = subfolder_list(contains({subfolder_list(:).name},'SW'));
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
folder_targets = vertcat(folder_cell{:});
%get the number of maps to do
map_num = size(folder_targets,1);
%% Load the gaussian models

%define the main file path
model_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';

%load the file
%define the file tags
file_tags = {'_rawClusters'};

%for all the tags
for tags = 1:length(file_tags)
    %load the most recent file from the overlap ones
    listing = dir(strcat(model_path,'\*',file_tags{tags},'.mat'));
    dates = datetime({listing.date});
    [~,ind] = max(dates);
    load_name = listing(ind).name;
    
    %load all the variables in the file
    load(fullfile(model_path,load_name))
end

%get some of the variables needed below
num_positions = max(trace2folder(:,1));
%number of cells
cell_num = max(trace2folder(:,4));
%number of maps
map_num = max(trace2folder(:,2));

%define the interval to grab. Each recording is a second and the desired
%range is from 100 to 300 ms
trace_range = 1000:3000;
%define the time interval for the background
trace_background = 1:999;
%define the target (synaptic) window
target_window = 71:1570;

%turn the throw away traces that were still over 3 stds into NaNs too
trace2folder(trace2folder(:,5)==0,5) = NaN;
%% Plot the clusters
close all
%for both polarities
for polars = 1:2
    %load the variables from the EI_cell
    clunum = EI_cell{8,polars};
    trace_ave = EI_cell{3,polars};
    trace_std = EI_cell{4,polars};
    clu_mem = EI_cell{14,polars};
    %generate a colormap for the cluster averages
    c_map = jet(clunum);
    %define the length of the trace to plot
    plot_length = 1:2001;
    figure
    %for all the clusters
    for clu = 1:clunum
%         figure
    %     plot(1:size(plot_matrix,1),trace_ave(clu,:)'+clu)
        %define the x vector
        x_vec = (0:size(plot_length,2)-1)./10;
        subplot(round(sqrt(clunum)),ceil(sqrt(clunum)),clu)
        shadedErrorBar(x_vec,trace_ave(clu,plot_length)'...
            ,trace_std(clu,plot_length)','lineprops',{'Color',c_map(clu,:)})
        hold('on')
        %plot the zero line
        plot(x_vec',zeros(length(x_vec),1),'k--')
        plot([7 7],get(gca,'YLim'),'k--')
        ylabel('Current (pA)')
        xlabel(strcat('Time (ms) #:',num2str(clu_mem(clu))))

    end
end
%% Generate maps based on the categorization from the clustering

%create a waitbar
w_bar = waitbar(0,'Calculating maps');

%allocate memory to store the cell's name and maps
all_maps = cell(cell_num,3);
%for all the cells
for cells = 1:cell_num
    %show progress
    waitbar(cells/cell_num,w_bar)
    
    %save the cell's name
    all_maps{cells,1} = uni_cells{cells};
    
    %get the info for the traces from this cell
    curr_cell = trace2folder(trace2folder(:,4)==cells,:);
    
    %for the polarities present
    for polars = 1:length(unique(curr_cell(:,3)))
        %get the set of maps present in this cell and polarity
        curr_maps = unique(curr_cell(curr_cell(:,3)==polars-1,2));
        %get the paths to the cells maps
        cell_paths = folder_all(curr_maps);
        %allocate memory to store the processed maps
        proc_maps = zeros(num_positions,length(cell_paths));
        %for all the maps
        for maps = 1:length(cell_paths)
            %fetch the background subtracted map
            bsub_map = trace_fetch(cell_paths{maps},trace_range,trace_background,num_positions);
            %also get the info for the current traces
            curr_traces = curr_cell(curr_cell(curr_cell(:,3)==polars-1&curr_cell(:,2)==curr_maps(maps)),:);
            %process the traces according to the type of response
            %for all the traces
            for trace = 1:num_positions
                %get the response type
                resp_type = curr_traces(trace,5);
                %if it's a NaN, skip the processing (i.e. blank trace)
                if isnan(resp_type)
                    continue
                else
                    %process the trace accordingly
                    proc_maps(trace,maps) = Trace_process(bsub_map(target_window,trace),resp_type,polars);
                end
            end
        end
        
        %blank the positions that only have a map in a single repetition
        %get the positions
%         single_blank = prod(proc_maps,2)==0;
%         proc_maps(single_blank,:) = 0;
        
        %save the average maps
        all_maps{cells,polars+1} = mean(reshape(proc_maps,sqrt(num_positions),sqrt(num_positions),[]),3);
        
    end
 
end
close(w_bar)
%% Plot finished maps
close all

figure
%initialize a counter for the subplots
plot_c = 1;
%select the target polarity
polars = 2;
%for all the cells
for cells = 1:cell_num
    
    %create a new figure every 10 cells
    if mod(cells,25)==0
        figure
        plot_c = 1;
    end
    
    subplot(5,5,plot_c)
    
    imagesc(all_maps{cells,polars+1})
    %update the plot counter
    plot_c = plot_c + 1;
    
end
%% Load full responses for plotting
%allocate memory to store the maps
temp_xsg = load(folder_all{1},'-mat');
%get the size of the maps
map_size = length(temp_xsg.data.ephys.trace_1);
%also the number of positions
num_positions = temp_xsg.header.mapper.mapper.positionNumber;
%define the interval to grab. Each recording is a second and the desired
%range is from 100 to 300 ms
trace_range = 1000:3000;
%define the time interval for the background
trace_background = 1:999;

%allocate a matrix for them
% map_matrix = zeros(length(trace_range),length(rand_ind));
map_matrix = zeros(length(trace_range)*num_positions,map_num);
%allocate memory for the standard deviation of the background for each
%trace
background_std = zeros(num_positions,map_num);

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
    %trim the trace
    temp_trace = temp_trace(trace_range,:);
    %subtract background activity
    temp_trace = temp_trace - background_act;
    %also load the map order to have all the maps ordered in their location
    %instead of in time of stimulation
    %linearize the map and apply to the data
    temp_trace = temp_trace(:,temp_map(:));
    %reshape again and store
    map_matrix(:,maps) = reshape(temp_trace,length(trace_range)*num_positions,[]);
end

close(w_bar)
%% Plot full maps with traces and selected response types


%% OFF Go map by map allocating the traces into one of the model categories

% %for all the maps
% for maps = 1:map_num
%     
%     %load the xsg file
%     temp_xsg = load(folder_targets{maps},'-mat');
%     %check the number of grid points in the cell
%     temp_map = temp_xsg.header.mapper.mapper.mapPatternArray;
%     %if the trace has less than 16x16
%     if size(temp_map,1) < 16 || size(temp_map,2) < 16 ||...
%             isempty(temp_xsg.data.ephys)
%         %mark the position and skip the trace
%         elim_vec(maps) = 0;
%         continue
%     end
%     %load the time trace
%     temp_trace = temp_xsg.data.ephys.trace_1;
%     %reshape to trim the relevant part of the data only
%     temp_trace = reshape(temp_trace,[],num_positions);
%    
%     %get the background activity
%     background_act = mean(temp_trace(trace_background,:),1);
%     %trim the trace
%     temp_trace = temp_trace(trace_range,:);
%     %subtract background activity
%     temp_trace = temp_trace - background_act;
%     %also load the map order to have all the maps ordered in their location
%     %instead of in time of stimulation
%     %linearize the map and apply to the data
%     temp_trace = temp_trace(:,temp_map(:));
%     %determine if it's exc or inh
%     map_polarity = (mean(temp_trace(:))>0) + 1;
%     
%     
%     
%     
% end
%% Save the maps excluding the direct responses