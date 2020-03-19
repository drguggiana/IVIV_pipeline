function PIPE_Trace_interpolate_joel(folder_all,ori_setup)

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
%     background_std(:,maps) = mean(std(temp_trace(trace_background,:),0,1));
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
        soma_center(maps,1) = soma_center(maps,1) - temp_xsg.header.mapper.mapper.xPatternOffset;
        soma_center(maps,2) = soma_center(maps,2) - temp_xsg.header.mapper.mapper.yPatternOffset;
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
excinh_path = 'R:\Share\Simon\Drago_Volker_Simon\Joel_Project\';

%load the file
load(fullfile(excinh_path,strcat('ExcInh_map_',ori_setup,'.mat')),'ExcInh_map');

%get only the binary values for the responses
map_polarity = vertcat(ExcInh_map{:,1});

% %calculate the average value for each map, since area under the curve
% %should report whether it is exc or inh
% ave_value = mean(map_matrix,1);

% figure
% histogram(ave_value)
% %calculate the polarity of each map to propagate later (0 is exc)
% map_polarity = ave_value>0;
%% Create maps for polarity and before/after

% %polarity map, 1 is excitation
% map_polarity = contains(folder_all,'inhibition');

%time map, 1 is pre (before)
time_map = contains(folder_all,'before');
%% Reshape the matrix, so that all traces are in a single column

map_matrix2 = reshape(map_matrix,length(trace_range),[]);
%linearize the std matrices also
background_std = background_std(:);

%and create a map from traces to folders (including polarity)
%allocate memory for the map
trace2folder = zeros(num_positions,size(map_matrix,2),6);
%and for the time_map separately
time_pertrace = zeros(num_positions,size(map_matrix,2));
%for all the folders
for folders = 1:size(folder_all,1)
    trace2folder(:,folders,1) = 1:num_positions;
    trace2folder(:,folders,2) = folders;
    trace2folder(:,folders,3) = map_polarity(folders);
    trace2folder(:,folders,4) = cell_id(folders);
    
    %fill in the time map too
    time_pertrace(:,folders) = time_map(folders);
    
end
%also reshape this map to be able to refer to the original map
trace2folder = reshape(trace2folder,[],6);
%also reshape the time map
time_pertrace = reshape(time_pertrace,[],1);
%% Filter out traces that are too flat (using 3 std criterion)
close all

%define the responsivity threshold
std_threshold = 3;
%define the windowing threshold
window_threshold = 3;


%calculate the std of each trace
std_matrix = std(map_matrix2,0,1);

%plot the distribution of the std
figure
histogram(std_matrix)

%get the direct and direct mix traces
direct_std = direct_std(:);
direct_map = direct_std>window_threshold.*background_std;

pseudodirect_std = pseudodirect_std(:);
pseudodirect_map = pseudodirect_std>window_threshold.*background_std;

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
prctile_map = std_matrix'>std_threshold.*background_std;
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
%% Bin the data

%temporally bin the data by a defined factor
%define the factor
bin_factor = 10;

%bin by the factor

%allocate memory for the binned data
bin_matrix = zeros(floor(size(map_matrix2,1)/bin_factor),size(map_matrix2,2));

%get the binning map
bin_map = discretize(1:size(map_matrix2,1),1:bin_factor:size(map_matrix2,1));

%for all the bins
for bins = 1:max(bin_map)
    %bin the data
    bin_matrix(bins,:) = mean(map_matrix2(bin_map==bins,:),1);
end
%% 3D interpolation to infer synaptic responses hidden within the direct responses

%create the waitbar
w_bar = waitbar(0,'Calculating maps');

%define a cell with the types of interpolation
interp_var = {'nearest'};
%allocate memory to store the interpolated maps
interp_cell = cell(map_num,length(interp_var));

%for all the maps
for maps = 1:map_num
    
    %update the waitbar
    waitbar(maps/map_num,w_bar)
    %get the info for the current map
    curr_info = trace2folder(trace2folder(:,2)==maps,:);
    %get the polarity
    polarity = mode(curr_info(:,3));
%     %if an inhibitory map, skip it
%     if polarity == 1
%         continue
%     end
    
    %check if there are any direct responses
    num_direct = sum(curr_info(:,5)==0);
    %if no direct responses, skip the interpolation and load the binned map
    if num_direct == 0
        interp_cell{maps,1} = bin_matrix(:,trace2folder(:,2)==maps);
        continue
    end
    
%     %don't interpolate the post TTX maps
%     if mode(time_pertrace(trace2folder(:,2)==maps))==0
%         continue
%     end
    
    %load the complete map
    curr_map = bin_matrix(:,trace2folder(:,2)==maps);

    %blank the direct traces
    curr_map(:,curr_info(:,5)==0) = NaN;
    
    %reshape the current map to a 3d grid (i.e. x , y and time)
    curr_map = reshape(curr_map,[],sqrt(num_positions),sqrt(num_positions));
    
%     figure
%     imagesc(squeeze(sum(curr_map,1)))
    
    %set up the interpolation
    
    %get the 3d coordinates of the direct responses
    direct_coord = find(isnan(curr_map));
    [qx,qy,qz] = ind2sub(size(curr_map),direct_coord);
    %and of the non-direct responses
    nondirect_coord = find(~isnan(curr_map));
    [dx,dy,dz] = ind2sub(size(curr_map),nondirect_coord);
    
    %for both kinds of interpolation
    for interpkind = 1:length(interp_var)
        filled_coord = griddatan([dx,dy,dz],curr_map(nondirect_coord),[qx,qy,qz],interp_var{interpkind});
        filled_map = curr_map;
        filled_map(direct_coord) = filled_coord;
        %store the interpolated map
        interp_cell{maps,interpkind} = reshape(filled_map,[],num_positions);
    end
%     figure
%     imagesc(squeeze(sum(filled_map,1)))
    
end
%close the waitbar
close(w_bar)
%% Replace the direct responses by 4 on the category column of trace2folder

trace2folder(trace2folder(:,5)==0,5) = 4;
%% Save the results

% %define the save path
% save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';
%define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out_joel';
%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_rawClass',ori_setup,'.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(out_path,save_name),'folder_all','prctile_map',...
    'trace_background','trace_range','uni_cells',...
    'trace2folder','background_std','trace_direct','trace_pseudodirect','soma_unique','interp_cell')
%terminate execution
return
%% Plot the maps

close all
%for all the maps
for maps = 1:map_num
    
    %if it's an empty cell, skip it
    if isempty(interp_cell{maps,1})
        continue
    end
    figure
    
    %get the current map
    curr_map = bin_matrix(:,trace2folder(:,2)==maps);
    %get the info associated with the map
    curr_info = trace2folder(trace2folder(:,2)==maps,:);
    
%     figure
%     imagesc(squeeze(sum(curr_map,1)))
%     title('Original image')
%     
%     figure
%     imagesc(squeeze(sum(filled_map,1)))
%     title('Filled image')
    

%     %turn curr map back into traces
%     curr_traces = reshape(curr_map,size(curr_map,1),16*16);
%     figure
%     plot(curr_traces)
    
    %compare the two maps
    %assemble a matrix with both maps
    plot_matrix = cat(4,permute(interp_cell{maps,1},[3 2 1]),permute(interp_cell{maps,2},[3 2 1]),...
        permute(reshape(curr_map,size(curr_map,1),16,16),[3 2 1]));
    %for both maps
    for submaps = 1:3
    
        plot_traces = plot_matrix(:,:,:,submaps);
        %define the amplitude factor
        amp_rat = 2;
        %define the subsampling factor
        sub_rat = 10;
        %define the separation factor
        sep_factor = (size(plot_traces,3) + 1000)/sub_rat;
        %for all the traces
        for x = 1:size(plot_traces,1)
            for y = 1:size(plot_traces,2)
                %get the index corresponding to this position in single index
                curr_ind_clu = sub2ind([size(plot_traces,1),size(plot_traces,2)],y,x);

                %otherwise make it blue
                switch submaps
                    case 1
                        trace_color = 'r';
                    case 2
                        trace_color = 'b';
                    case 3
                        if curr_info(sub2ind([16,16],y,x),5)==0
                            trace_color = 'k';
                        else
                            trace_color = 'm';
                        end
                end
                %                     end
                %get the x vector, correcting for the array position
                x_vec = (1:length(squeeze(plot_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
                %and the y vector
                y_vec = squeeze(plot_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
                %plot the result
                plot(x_vec,-y_vec,trace_color)
                hold('on')
                %plot a 0 line
                plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
                %plot a line at 7 ms
                plot(x_vec([7 7]),-[max(y_vec) min(y_vec)],'g--')

            end
        end
    end
end
%% Average the reps and compare pre and post TTX

%allocate memory to store the averaged maps (post exc, pre ext, interp exc
%lin and nearest ,then same for inh in last dimension)
map_cell = cell(max(trace2folder(:,4)),4,2);
%allocate memory to store the info from the averaged maps
map_info = cell(max(trace2folder(:,4)),2);
%for all the cells
for cells = 1:max(trace2folder(:,4))
    
    %for both polarities
    for polarity = 0:1
        %for before and after
        for timet = 0:1
            %get the idx for each of the for loop iterations
            cell_idx = trace2folder(:,4)==cells;
            polar_idx = trace2folder(:,3)==polarity;
            time_idx = time_pertrace==timet;
            %get the info for these conditions
            curr_info = trace2folder(cell_idx&polar_idx&time_idx,:);
            %get the id of the current maps
            curr_ids = unique(curr_info(:,2));
            %get the number of maps in the cells
            curr_num = length(curr_ids);
            %allocate memory to store the average map
            curr_map = zeros(size(bin_matrix,1),num_positions);
            %allocate temporary memory to store the trace info
            temp_info = zeros(num_positions,size(trace2folder,2),curr_num);
            %for all the maps
            for submaps = 1:curr_num
                %get the idx of the current map
                map_idx = trace2folder(:,2)==curr_ids(submaps);
                %get the current map and add it for averaging
                curr_map = curr_map + bin_matrix(:,cell_idx&polar_idx&time_idx&map_idx)./curr_num;
                %store the trace info
                temp_info(:,:,submaps) = trace2folder(cell_idx&polar_idx&time_idx&map_idx,:);
            end
            %store the average in the corresponding cell
            map_cell{cells,timet+1,polarity+1} = curr_map;
            %if a before map
            if timet == 1
                %save the info for the map
                map_info{cells,polarity+1} = mode(temp_info,3);
            end
            
            
        end
        %do the same for the interpolated map
        %allocate memory to store the average map
        interp_map = zeros(size(bin_matrix,1),num_positions,2);
        
        %for all the maps
        for submaps = 1:curr_num
            %get the idx of the current map
            map_idx = trace2folder(:,2)==curr_ids(submaps);
           
            %for both interpolation methods
            for intermethod = 1:2
                %if the slot is empty, skip it
                if isempty(interp_cell{curr_ids(submaps),intermethod})
                    continue
                end
                %get the current map and add it for averaging
                interp_map(:,:,intermethod) = nansum(cat(4,interp_map(:,:,intermethod),...
                    reshape(interp_cell{curr_ids(submaps),intermethod},[],num_positions)./curr_num),4);
            end
        end
        
        %store in the corresponding cell
        map_cell{cells,3,polarity+1} = interp_map(:,:,1);
        map_cell{cells,4,polarity+1} = interp_map(:,:,2);
        
        
    end
    
end
%% Plot the overlap maps

close all

%for all the maps
for cells = 1:max(trace2folder(:,4))
    
    %for both polarities
    for polarity = 1
        figure
        
%         %get the current map
%         curr_map = bin_matrix(:,trace2folder(:,2)==maps);
%         %get the info associated with the map
%         curr_info = trace2folder(trace2folder(:,2)==maps,:);
        
        
        %compare the two maps
        %assemble a matrix with both maps
        plot_matrix = cat(3,map_cell{cells,:,polarity});
        %for both maps
        for submaps = 1:4
            
            plot_traces = permute(reshape(plot_matrix(:,:,submaps),[],sqrt(num_positions),sqrt(num_positions)),[3 2 1]);
            %define the amplitude factor
            amp_rat = 2;
            %define the subsampling factor
            sub_rat = 10;
            %define the separation factor
            sep_factor = (size(plot_traces,3) + 1000)/sub_rat;
            %for all the traces
            for x = 1:size(plot_traces,1)
                for y = 1:size(plot_traces,2)
                    %get the index corresponding to this position in single index
                    curr_ind_clu = sub2ind([size(plot_traces,1),size(plot_traces,2)],y,x);
                    
                    %otherwise make it blue
                    switch submaps
                        case 1
                            trace_color = 'r';
                        case 2
                            trace_color = 'b';
                        case 3
                            trace_color = 'k';
                        case 4
                            trace_color = 'm';
                            
                    end
                    %                     end
                    %get the x vector, correcting for the array position
                    x_vec = (1:length(squeeze(plot_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
                    %and the y vector
                    y_vec = squeeze(plot_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
                    %plot the result
                    plot(x_vec,-y_vec,trace_color)
                    hold('on')
                    %plot a 0 line
                    plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
                    %plot a line at 7 ms
                    plot(x_vec([7 7]),-[max(y_vec) min(y_vec)],'g--')
                    
                end
            end
        end
    end
end
%% Attempt to recreate mixed responses by convolving interpolated synaptic and post TTX

%allocate memory for the convolved traces (for both polarities)
conv_cell = cell(max(trace2folder(:,4)),2);
%for all the maps
for cells = 1:max(trace2folder(:,4))
    %for both polarities
    for polarity = 1
        %load the trace info
        trace_info = map_info{cells,polarity};
        %load the corresponding traces
        post_traces = map_cell{cells,1,polarity};
        nearest_traces = map_cell{cells,4,polarity};
        pre_traces = map_cell{cells,2,polarity};
        %allocate memory for the output
        conv_traces = zeros(size(post_traces));
        %for all the traces
        for traces = 1:num_positions
            %if the trace is not a direct response, skip it
            if trace_info(traces,5)~= 0
                %copy the pre trace in this position
                conv_traces(:,traces) = pre_traces(:,traces);
                continue
            end
           
            %convolve the post and nearest traces to generate a mixed
            %response
%             conv_traces(:,traces) = normr_2(conv(nearest_traces(:,traces),post_traces(:,traces),'same')).*min((pre_traces(:,traces)));
%             conv_traces(:,traces) = -(normr_2(sum([nearest_traces(:,traces),post_traces(:,traces)],2)).*min(pre_traces(:,traces)))+min(pre_traces(:,traces));
            conv_traces(:,traces) = sum([nearest_traces(:,traces),post_traces(:,traces)],2);
        end
        %store in the main convolution cell
        conv_cell{cells,polarity} = conv_traces;
    end
end
%% Plot the convolution results
close all

%for all the maps
for cells = 1%:max(trace2folder(:,4))
    
    %for both polarities
    for polarity = 1
        figure
        
%         %get the current map
%         curr_map = bin_matrix(:,trace2folder(:,2)==maps);
%         %get the info associated with the map
%         curr_info = trace2folder(trace2folder(:,2)==maps,:);
        
        
        %compare the two maps
        %assemble a matrix with both maps (first original, then
        %convolution)
        plot_matrix = cat(3,map_cell{cells,2,polarity},conv_cell{cells,polarity},map_cell{cells,1,polarity});
        %for both maps
        for submaps = 1:3
            
            plot_traces = permute(reshape(plot_matrix(:,:,submaps),[],sqrt(num_positions),sqrt(num_positions)),[3 2 1]);
            %define the amplitude factor
            amp_rat = 2;
            %define the subsampling factor
            sub_rat = 10;
            %define the separation factor
            sep_factor = (size(plot_traces,3) + 1000)/sub_rat;
            %for all the traces
            for x = 1:size(plot_traces,1)
                for y = 1:size(plot_traces,2)
                    %get the index corresponding to this position in single index
                    curr_ind_clu = sub2ind([size(plot_traces,1),size(plot_traces,2)],y,x);
                    
                    %otherwise make it blue
                    switch submaps
                        case 1
                            trace_color = 'r.-';
                        case 2
                            trace_color = 'b.-';
                        case 3
                            trace_color = 'k';
%                         case 4
%                             trace_color = 'm';
                    end
                    %                     end
                    %get the x vector, correcting for the array position
                    x_vec = (1:length(squeeze(plot_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
                    %and the y vector
                    y_vec = squeeze(plot_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
                    %plot the result
                    plot(x_vec,-y_vec,trace_color)
                    hold('on')
                    %plot a 0 line
                    plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
                    %plot a line at 7 ms
                    plot(x_vec([7 7]),-[max(y_vec) min(y_vec)],'g--')
                    
                end
            end
        end
    end
end