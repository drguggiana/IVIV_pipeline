%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load paths

%pick the folders to use
folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\INPUT MAPS_final\Setup2_TTX');


%define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';
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
folder_all = vertcat(folder_cell{:});
%get the number of maps to do
map_num = size(folder_all,1);
%% Create maps for polarity and before/after

%polarity map, 1 is excitation
polar_map = contains(folder_all,'excitation');

%time map, 1 is before
time_map = contains(folder_all,'before');
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
%get the size of the maps
map_size = length(temp_xsg.data.ephys.trace_1);
%also the number of positions
num_positions = temp_xsg.header.mapper.mapper.positionNumber;
%define the interval to grab. Each recording is a second and the desired
%range is from 100 to 300 ms
trace_range = 1001:3000;
%define the time interval for the background
trace_background = 1:1000;

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

%get rid of the empty spaces with the skipped traces
map_matrix = map_matrix(:,elim_vec==1);
%also within the std matrix
background_std = background_std(:,elim_vec==1);
% map_matrix2 = map_matrix;
%and also modify the folder vector (so I can refer back to the particular
%maps)
folder_all = folder_all(elim_vec==1);
%update the cell name array
cell_names = cell_names(elim_vec==1);

%get the unique cell names
[uni_cells,~,cell_id] = unique(cell_names);
%% Reshape the matrix, so that all traces are in a single column

map_matrix2 = reshape(map_matrix,length(trace_range),[]);
%linearize the std matrix also
background_std = background_std(:);
%and create a map from traces to folders (including polarity)
%allocate memory for the map
trace2folder = zeros(num_positions,size(map_matrix,2),6);
%for all the folders
for folders = 1:size(folder_all,1)
    trace2folder(:,folders,1) = 1:num_positions;
    trace2folder(:,folders,2) = folders;
    trace2folder(:,folders,3) = polar_map(folders);
    trace2folder(:,folders,4) = time_map(folders);
    trace2folder(:,folders,5) = cell_id(folders);
end
%also reshape this map to be able to refer to the original map
trace2folder = reshape(trace2folder,[],6);
%% Filter out traces that are too flat (using 3 std criterion)
close all
%calculate the std of each trace
std_matrix = std(map_matrix2,0,1);
% %define the percentile threshold
% prctile_thres = 80;
% prctile_cutoff = prctile(std_matrix,prctile_thres);
%plot the distribution of the std
figure
histogram(std_matrix)
% %exclude traces on the lowest 10th percentile
% map_matrix2 = map_matrix2(:,std_matrix>prctile_cutoff);
% %exclude traces with less than 3 stds over background
% map_matrix2 = map_matrix2(:,std_matrix'>3.*background_std);
%also save a map of the selected traces
prctile_map = std_matrix'>3.*background_std;
%generate a map of the origin of the traces left after the percentile cut
postperc_map = trace2folder(prctile_map,:);
%mark the places on the trace2folder matrix
trace2folder(prctile_map==0,6) = NaN;
%% Compare traces from before and after

close all
%get the trace length
trace_length = length(map_matrix(:,1))/num_positions;
%allocate memory to store the subtraction
sub_cell = cell(2,1);
%for all the cells
for cells = 1:max(trace2folder(:,5))
    %define polarity
    polarity = 0;
    
    switch polarity
        case 0
            pol_str = 'Inh';
        case 1
            pol_str = 'Exc';
    end
    figure
    %for before, after and subtraction
    for time_c = 0:1
        
        switch time_c
            case 0
                time_str = 'After';
            case 1
                time_str = 'Before';
        end
        
        switch time_c
            case {0,1}
                %get the target map
                target_map = unique(trace2folder(trace2folder(:,3)==polarity&trace2folder(:,5)==cells&trace2folder(:,4)==time_c,2));
                %allocate memory to store the averaged maps
                map_traces = zeros(sqrt(num_positions),sqrt(num_positions),trace_length);
                %for all the maps
                for maps = target_map'
                    
                    %get the map from the main matrix
                    map_traces = map_traces + permute(reshape(map_matrix(:,maps),...
                        [],sqrt(num_positions),sqrt(num_positions)),[3 2 1])/length(target_map);
                    
                end
                %store the result for subtracting later
                sub_cell{time_c+1} = map_traces;
            case 2
                %calculate the subtraction of the traces
                map_traces = sub_cell{2}-sub_cell{1};
                
                
        end
        
        %define the amplitude factor
        amp_rat = 2;
        %define the subsampling factor
        sub_rat = 10;
        %define the separation factor
        sep_factor = (size(map_traces,3) + 1000)/sub_rat;
        %for all the traces
        for x = 1:size(map_traces,1)
            for y = 1:size(map_traces,2)
                %get the index corresponding to this position in single index
                curr_ind_clu = sub2ind([size(map_traces,1),size(map_traces,2)],y,x);
   
                %otherwise make it blue
                switch time_c
                    case 0
                        trace_color = 'r';
                    case 1
                        trace_color = 'b';
                    case 2
                        trace_color = 'm';
                end
                %                     end
                %get the x vector, correcting for the array position
                x_vec = (1:length(squeeze(map_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
                %and the y vector
                y_vec = squeeze(map_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
                %plot the result
                plot(x_vec,-y_vec,trace_color)
                hold('on')
                %plot a 0 line
                plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
                %plot a line at 7 ms
                plot(x_vec([7 7]),-[max(y_vec) min(y_vec)],'g--')
               
            end
        end
        title(strcat('Cell: ',uni_cells{cells},'Pol:',pol_str))
         
    end
end
%% Absolute error calculation

close all

%get the delta between before and after

cell_num = max(trace2folder(:,5));
%allocate memory to store the deltas (for each polarity and delta vs amplitude)
delta_cell = cell(cell_num,2,2);
%also allocate memory to keep the paired traces (i.e. can be averaged in
%both conditions and are guided by the before)
paired_maps = cell(cell_num,2);

%for all the cells
for cells = 1:cell_num
    %for both polarities
    for polarity = 0:1
        %find if the cell has before and after for this polarity, otherwise
        %skip the loop iteration
        if length(unique(trace2folder(trace2folder(:,5)==cells&trace2folder(:,3)==polarity,4))) < 2
            continue
        end
        %get the different sets of indexes
        cell_idx = trace2folder(:,5)==cells;
        polarity_idx = trace2folder(:,3)==polarity;
        active_idx = ~isnan(trace2folder(:,6));
  
        %allocate memory to store the xor positions
        xor_cell = cell(2,1);
        %for both times
        for times = 0:1
            %get the positions that have active traces in the condition
            positions = trace2folder(cell_idx&polarity_idx&active_idx&trace2folder(:,4)==times,1);

            %get the map of origin of each position
            ori_map = trace2folder(cell_idx&polarity_idx&active_idx&trace2folder(:,4)==times,2);

            %get the positions that are present in all repeats
            %get the unique folders in the set
            rep_vec = unique(ori_map);
            %get the number of reps
            rep_num = length(rep_vec);
            %allocate memory to store the rep positions
            rep_cell = cell(rep_num,1);
            for reps = 1:rep_num
                rep_cell{reps} = positions(ori_map==rep_vec(reps));
            end
            %if two reps
            if rep_num == 2
                %get the intersection between the reps
                xor_cell{times+1} = setxor(rep_cell{1},rep_cell{2});
            %if 1 rep, just save the cell itself
            elseif rep_num == 1 
                xor_cell{times+1} = rep_cell{1};
            end
        end
        %get the unique positions to blank
        xor_pos = unique(vertcat(xor_cell{:}));
        %get the positions that have active traces in the before condition,
        %to use them in the after condition also
        positions = trace2folder(cell_idx&polarity_idx&active_idx&trace2folder(:,4)==1,1);
        %exclude the positions that can't be averaged
        common_pos = setxor(positions,xor_pos);
        %get the number of positions
        common_num = length(common_pos);
        %calculate the average of all the positions
        %allocate memory for the average and for the background std
        ave_pos = cell(2,3);
        %for all the times
        for times = 0:1
            %get the traces
            curr_traces = map_matrix2(:,cell_idx&polarity_idx&trace2folder(:,4)==times);
            %and the background std
            curr_bckg = background_std(cell_idx&polarity_idx&trace2folder(:,4)==times);
            %and the current activity
            curr_act = trace2folder(cell_idx&polarity_idx&trace2folder(:,4)==times,6);
            %allocate memory for the average matrix for these times
            ave_mat = zeros(size(curr_traces,1),common_num);
            std_mat = zeros(common_num,1);
            active_mat = zeros(common_num,1);
            %for all the common positions
            for pos = 1:length(common_pos)
                %get the positions vector
                pos_vector = trace2folder(cell_idx&polarity_idx&trace2folder(:,4)==times,1)==common_pos(pos);
                %calculate the average
                ave_mat(:,pos) = mean(curr_traces(:,pos_vector),2);
                %also of the background std
                std_mat(pos) = mean(curr_bckg(pos_vector),1);
                %and get the activity of the position
                active_mat(pos) = any(~isnan(curr_act(pos_vector)));
            end
            %store the results in the storage cell
            ave_pos{times+1,1} = ave_mat;
            ave_pos{times+1,2} = std_mat;
            ave_pos{times+1,3} = active_mat;
        end
        
        %save the paired maps for later use
        paired_maps{cells,polarity+1} = ave_pos;
      
        %calculate the deltas between the corresponding positions
        delta_cell{cells,polarity+1,1} = mean(diff(cat(3,ave_pos{:,1}),1,3).^2,1);
        delta_cell{cells,polarity+1,2} = sum(ave_pos{2,1},1);
    end
end

%plot the results
%define the polarity labels
pol_label = {'Inh','Exc'};
%for each polarity
for polarity = 0:1
    figure
    %concatenate the cells
    all_diff = horzcat(delta_cell{:,polarity+1,1});
    all_amp = horzcat(delta_cell{:,polarity+1,2});
    
    
    
    histogram(all_diff)
%     set(gca,'XScale','log')
    title(pol_label{polarity+1})
    xlabel('Average square difference between Before and After')
    %plot the relationship between amplitude and error
    figure
    plot(all_amp,all_diff,'*')
    xlabel('Current')
    ylabel('Error')
    title(pol_label{polarity+1})
%     set(gca,'XScale','log')


end
%% Calculate the number of states for each map

%allocate memory to store the map variances (cell, polarity, times)
map_variance = cell(cell_num,polarity,2);

%for all the cells
for cells = 1:cell_num
    %for both polarities
    for polarity = 1:2
        %get the current maps
        curr_maps = paired_maps{cells,polarity};
        %skip the cell if it's empty
        if isempty(curr_maps)
            continue
        end
        
        %for both times
        for times = 1:2
            map_variance{cells,polarity,times} = std(curr_maps{times,1},0,1);
        end
        
    end
end
%% Sliding threshold contamination calculation

%create a matrix with all the active before positions and their after
%counterparts

%using the matrix, calculate now the proportion of direct vs synaptic
%traces per cell and also the amount of direct response incurred for a
%sliding window

%define the sliding window (in ms*10)
sliding_window = 1:1:100;

%get the number of positions
window_num = length(sliding_window);

%allocate memory to store the number of direct responses and the percentage
%variance (window,cell,polarity,count/current)
window_cell = cell(window_num,cell_num,2,3);

%for all the cells
for cells = 1:cell_num
    %for both polarities
    for polarity = 1:2
        %get the current maps
        curr_maps = paired_maps{cells,polarity};
        %skip the cell if it's empty
        if isempty(curr_maps)
            %for all the positions
            for winds = 1:window_num
                window_cell{winds,cells,polarity,1} = NaN;
                window_cell{winds,cells,polarity,2} = NaN;
                window_cell{winds,cells,polarity,3} = NaN;
            end
            continue
        end
        
%         %quantify the map variance in the before condition
%         before_current = sum(curr_maps{1,1});
%         %get the number of active positions in the after condition
%         active_pos = sum(curr_maps{1,3});
        %get the number of active positions in the before condition
        total_pos = sum(curr_maps{2,3});
        %calculate the total variance per map
        total_var = sum(map_variance{cells,polarity,2});
        %for all the positions
        for winds = 1:window_num
            %filter the traces by the window
            %define the current window
            curr_wind = sliding_window(winds);
            %calculate the std of the window
            window_std = std(curr_maps{2,1}(1:curr_wind,:),0,1)';
            %get the trace idx that get included based on the threshold
            filtered_traces = ~(window_std>3.*curr_maps{2,2});
            %store the percentage of after traces that pass the threshold
%             window_cell{winds,cells,polarity,1} = sum(curr_maps{1,3}&filtered_traces)/active_pos;
%             %calculate purity index: out of the total traces, what proportion has an
%             %active direct response
%             window_cell{winds,cells,polarity,1} = (sum(curr_maps{2,3}&filtered_traces)-sum(curr_maps{1,3}...
%                 &filtered_traces))/sum(curr_maps{2,3}&filtered_traces);
            %calculate purity index: out of the total traces, what proportion has an
            %active direct response
%             window_cell{winds,cells,polarity,1} = (sum(curr_maps{1,3})-sum(curr_maps{1,3}...
%                 &filtered_traces))/sum(curr_maps{1,3});
            window_cell{winds,cells,polarity,1} = sum(curr_maps{1,3}...
                &~filtered_traces)/sum(curr_maps{1,3});

            %also calculate the amount of the total current calculated due
            %to the direct responses
            window_cell{winds,cells,polarity,2} = sum(sum(abs(curr_maps{1,1}(curr_wind+1:end,filtered_traces))))/...
                sum(sum(abs(curr_maps{2,1}(curr_wind+1:end,filtered_traces))));
            %calculate the proportion of total variance kept by the traces
            window_cell{winds,cells,polarity,3} = sum(map_variance{cells,polarity,2}(filtered_traces))/total_var;
        end
    end
end



%plot the curves
close all
%get a set of colors for the cells
c_map = parula(cell_num);

%turn the window cell into an array
window_mat = cell2mat(window_cell);
%calculate average and std across cells
cell_ave = squeeze(nanmean(window_mat,2));
cell_std = squeeze(nanstd(window_mat,0,2));
%for both polarities
for polarity = 1:2
    switch polarity
        case 1
            pol_str = 'Inh';
        case 2
            pol_str = 'Exc';
    end
    figure('Name',pol_str)
    
    %for all the cells
    for cells = 1:cell_num
        %if the cell is empty, skip it
        if isempty(cat(1,window_cell{:,cells,polarity,:}))
            continue
        end
        %concatenate along the window dimensions
%         subplot(1,2,1)
        yyaxis left
        plot(sliding_window/10,((window_mat(:,cells,polarity,1))),'-','Color',c_map(cells,:))
        hold('on')
%         ylabel('Proportion of direct traces added')
        ylabel('Proportion of direct traces removed')
        yyaxis right
        plot(sliding_window(2:end)/10,cumsum((diff(window_mat(:,cells,polarity,3)))/sum(diff(window_mat(:,cells,polarity,3)))),'--','Color',c_map(cells,:))
        hold('on')
        ylabel('Proportion of variance covered')
        xlabel('Window size (ms)')
%         subplot(1,2,2)
%         plot(sliding_window/10,window_mat(:,cells,polarity,2))
%         hold('on')
%         xlabel('Window size (ms)')
%         ylabel('Proportion of direct signal present')
    end
    
%     %also plot the averages
%     subplot(1,2,2)
    yyaxis left
% errorbar(sliding_window/10,cell_ave(:,polarity,1),cell_std(:,polarity,1))
plot(sliding_window/10,cell_ave(:,polarity,1),'*')

%     shadedErrorBar(sliding_window/10,cell_ave(:,polarity,1),cell_std(:,polarity,1))
    yyaxis right
%     shadedErrorBar(sliding_window/10,cell_ave(:,polarity,3),cell_std(:,polarity,3))
% errorbar(sliding_window/10,cell_ave(:,polarity,3),cell_std(:,polarity,3))
plot(sliding_window/10,cell_ave(:,polarity,3),'*')
end