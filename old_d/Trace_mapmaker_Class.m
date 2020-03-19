%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load the pre-processed traces

%define the main file path
model_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';

%load the file
%define the file tags
file_tags = {'_rawClassSetup2'};

%for all the tags
for tags = 1:length(file_tags)
    %load the most recent file from the overlap ones
    listing = dir(strcat(model_path,'\*',file_tags{tags},'*.mat'));
    dates = datetime({listing.date});
    [~,ind] = max(dates);
    load_name = listing(ind).name;
    
    %load all the variables in the file
    load(fullfile(model_path,load_name))
end
%get the setup of origin
ori_setup = file_tags{1}(end-5:end);
%get some of the variables needed below
num_positions = max(trace2folder(:,1));
%number of cells
cell_num = max(trace2folder(:,4));
%number of maps
map_num = max(trace2folder(:,2));

% %define the interval to grab. Each recording is a second and the desired
% %range is from 100 to 300 ms
% trace_range = 1001:3000;
% %define the time interval for the background
% trace_background = 1:1000;


% %turn the throw away traces that were still over 3 stds into NaNs too
% trace2folder(trace2folder(:,5)==0,5) = NaN;
%% Generate maps based on the categorization from the clustering

%create a waitbar
w_bar = waitbar(0,'Calculating maps');

%fields on the info matrix
%EXC:1 peak direct 2 #direct, 3 #second window, 4 #synaptic
%INH:1 %pos/total 2 #direct and neg bigger, 3 #direct, pos bigger, 4
%#synaptic neg bigger,5 %  synaptic, pos bigger
%allocate memory to store the cell's name and maps
all_maps = cell(cell_num,8);
%for all the cells
for cells = 1:cell_num
    %show progress
    waitbar(cells/cell_num,w_bar)
    
    %save the cell's name
    all_maps{cells,1} = uni_cells{cells};
    %also the cells soma coordinates
    all_maps{cells,8} = soma_unique(cells,:);
    
    %get the info for the traces from this cell
    curr_cell = trace2folder(trace2folder(:,4)==cells,:);
    
    %for the polarities present
    for polars = 1:2
        

            %get the set of maps present in this cell and polarity
            curr_maps = unique(curr_cell(curr_cell(:,3)==polars-1,2));
            %get the paths to the cells maps
            cell_paths = folder_all(curr_maps);
            %allocate memory to store the processed maps
            proc_maps = zeros(num_positions,length(cell_paths));
            %allocate memory to store the 16x16x6 trace2folder array
            pixel_sub = zeros(num_positions,size(trace2folder,2),length(cell_paths));

            %for all the maps
            for maps = 1:length(cell_paths)
                %also get the info for the current traces
                curr_traces = curr_cell(curr_cell(:,3)==polars-1&curr_cell(:,2)==curr_maps(maps),:);
                %if there are interpolated traces, load the corresponding
                %map
                if exist('interp_cell','var')
                    %define the target (synaptic) window
                    target_window = 8:157;
                    %load the map from the interpolation cell
                    bsub_map = interp_cell{curr_maps(maps)};
                else
                    %define the target (synaptic) window
                    target_window = 71:1570;
                    %fetch the background subtracted map from the raw data
                    bsub_map = trace_fetch(cell_paths{maps},trace_range,trace_background,num_positions);
                end
                %count the types of responses for this map
                map_info = response_counter(curr_traces,polars,bsub_map);
                %store the trace2folder info
                pixel_sub(:,:,maps) = curr_traces;
                
                %process the traces according to the type of response
                %for all the traces
                for trace = 1:num_positions
                    %get the response type
                    resp_type = curr_traces(trace,5);
                    %if it's a NaN, skip the processing (i.e. blank trace with a zero)
                    if isnan(resp_type)
                        continue
                    end

                    %process the trace accordingly
                    proc_maps(trace,maps) = Trace_process(bsub_map(target_window,trace),resp_type,polars);
                    
                end                
                
            end

%             %blank the positions that only have a map in a single repetition
%             %get the positions
%             single_blank = prod(proc_maps,2)==0;
%             proc_maps(single_blank,:) = 0;

            %save the average maps
            all_maps{cells,polars+1} = nanmean(reshape(proc_maps,sqrt(num_positions),sqrt(num_positions),[]),3);
            %store the extra info
            all_maps{cells,polars+3} = map_info;
            %and also the pixel by pixel info in an array
            all_maps{cells,polars+5} = reshape(mode(pixel_sub,3),sqrt(num_positions),sqrt(num_positions),[]);
    end
 
end
close(w_bar)
%% Plot finished maps
close all

% %select the target polarity
% polars = 1;

%for both polarities
for polars = 1:2
    
    %initialize a counter for the subplots
    plot_c = 1;
    %for all the cells
    for cells = 1:cell_num

        %create a new figure every 10 cells
        if mod(cells,25)==1
            figure
            plot_c = 1;
        end

        subplot(5,5,plot_c)

        imagesc(all_maps{cells,polars+1},'AlphaData',~isnan(all_maps{cells,polars+1}))
        set(gca,'XTick',[],'YTick',[])
        %update the plot counter
        plot_c = plot_c + 1;

    end
end
%% Plot the composition of the maps
% close all

%for polarity
for polars = 1:2
    %get the combined matrices from all maps
    combo_mats = cat(2,all_maps{:,polars+3});
    figure
    if polars == 1
        plot_num = 4;
        plot_titles = {'Maximum deflection (pA)','# direct responses','# 2nd window responses','# synaptic responses'};
    else
        plot_num = 5;
        plot_titles = {'Percentage positive charge','# direct neg','# synaptic pos','# synaptic neg','# direct pos'};
    end
    %for all the plots
    for plots = 1:plot_num
        
        subplot(ceil(sqrt(plot_num)),round(sqrt(plot_num)),plots)
        histogram(combo_mats(plots,:))
        title(plot_titles{plots})
    end
    
end
%% OFF Load full responses for plotting
% %allocate memory to store the maps
% temp_xsg = load(folder_all{1},'-mat');
% %get the size of the maps
% map_size = length(temp_xsg.data.ephys.trace_1);
% %also the number of positions
% num_positions = temp_xsg.header.mapper.mapper.positionNumber;
% %define the interval to grab. Each recording is a second and the desired
% %range is from 100 to 300 ms
% trace_range = 1001:3000;
% %define the time interval for the background
% trace_background = 1:1000;
% 
% %allocate a matrix for them
% % map_matrix = zeros(length(trace_range),length(rand_ind));
% map_matrix = zeros(length(trace_range)*num_positions,map_num);
% % %allocate memory for the standard deviation of the background for each
% % %trace
% % background_std = zeros(num_positions,map_num);
% 
% %set up a vector for the skipped files
% elim_vec = ones(size(map_matrix,2),1);
% 
% %create a waitbar
% w_bar = waitbar(0,'Loading maps');
% %load the random maps into the matrix
% % for maps = 1:length(rand_ind)
% for maps = 1:map_num
%     
%     %show progress
%     waitbar(maps/map_num,w_bar)
%     %load the xsg file
% %     temp_xsg = load(folder_all{rand_ind(maps)},'-mat');
%     temp_xsg = load(folder_all{maps},'-mat');
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
% %     %get the std deviation of the background
% %     background_std(:,maps) = std(temp_trace(trace_background,:),0,1);
%     %trim the trace
%     temp_trace = temp_trace(trace_range,:);
%     %subtract background activity
%     temp_trace = temp_trace - background_act;
%     %also load the map order to have all the maps ordered in their location
%     %instead of in time of stimulation
%     %linearize the map and apply to the data
%     temp_trace = temp_trace(:,temp_map(:));
%     %reshape again and store
%     map_matrix(:,maps) = reshape(temp_trace,length(trace_range)*num_positions,[]);
% end
% 
% close(w_bar)
%% OFF Plot full maps with traces and selected response types

% close all
% 
% %define the cluster to focus on
% resp_type = [1 2 3 4];
% %set whether to look at exc (0) or inh (1)
% exc_inh = 0;
% 
% % %list of dirty exc maps
% % dirty_exc = [23 26 27 62 93 98 137 138 155 160 173];
% % %set whether to look at synaptic or direct responses
% % character = 1;
% % %load the corresponding cluster indexes
% % clusters = EI_cell{13,exc_inh+1};
% %get the IDs of the traces corresponding to this cluster
% % trace_id = temp_postperc_map(clusters==clu_focus,:);
% % trace_id = trace2folder(any(trace2folder(:,5)==resp_type,2)&trace2folder(:,3)==exc_inh,:);
% trace_id = trace2folder(trace2folder(:,3)==exc_inh,:);
% % %also get the corresponding cluster indexes
% % temp_idx = trace2folder(trace2folder(:,3)==exc_inh&~isnan(trace2folder(:,5)),5);
% % clu_idx = clusters(any(temp_idx==resp_type,2));
% %get the map ids for all traces in this cluster
% map_ids = unique(trace_id(:,2));
% 
% %for all the maps
% % for target_map = map_ids(dirty_exc)'%
% for target_map = map_ids([20 21])'%map_ids(randperm(length(map_ids),20))'
%     %define the map of interest
% %     target_map = 60;
%     figure
%     %get the map from the main matrix
%     map_traces = permute(reshape(map_matrix(:,target_map),[],sqrt(num_positions),sqrt(num_positions)),[3 2 1]);
% 
%     %define the amplitude factor
%     amp_rat = 2;
%     %define the subsampling factor
%     sub_rat = 10;
%     %define the separation factor
%     sep_factor = (size(map_traces,3) + 1000)/sub_rat;
%     %for all the traces
%     for x = 1:size(map_traces,1)
%         for y = 1:size(map_traces,2)
%             %get the index corresponding to this position in single index
% %             curr_ind = sub2ind([size(map_traces,1),size(map_traces,2)],x,y);
%             curr_ind_clu = sub2ind([size(map_traces,1),size(map_traces,2)],y,x);
%             
%             %select the trace color depending on the response type
%             trace_resp = trace2folder(trace2folder(:,1)==curr_ind_clu & trace2folder(:,2)==target_map,5);
%             switch trace_resp
%                 case 0
%                     trace_color = 'r';
%                 case 1
%                     trace_color = 'b';
%                 case 2
%                     trace_color = 'm';
%                 otherwise
%                     trace_color = 'k';
%             end
% %             %if the trace belongs to the target color
% %             if any(trace_id(:,1)==curr_ind_clu & trace_id(:,2)==target_map)
% %                 %make it red
% %                 trace_color = 'r';
% %             else
% %                 %otherwise make it blue
% %                 trace_color = 'b';
% %             end
%             %get the x vector, correcting for the array position
%             x_vec = (1:length(squeeze(map_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
%             %and the y vector
%             y_vec = squeeze(map_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
%             %plot the result
%             plot(x_vec,-y_vec,trace_color)
%             hold('on')
%             %plot a 0 line
%             plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
%             %plot a line at 7 ms
%             plot(x_vec([7 7]),-[max(y_vec) min(y_vec)],'g--')
% %             %plot the cluster of origin
% %             text(x_vec(1),-y_vec(1),num2str(clu_idx(trace_id(:,1)==curr_ind_clu &trace_id(:,2)==target_map)))
%             
%         end
%     end
% end
%% OFF Given a list of maps, get the cell names

% %define the list of maps (in this case, exc ones with positive deflections)
% map_list = [45 48 49 119 176 184 256 257 288 296 319];
% %get the number of maps
% list_num = length(map_list);
% %allocate memory for the list of cell names
% name_list = cell(list_num,2);
% %for all the maps on the list
% for maps = 1:list_num
%     %get the id number of the cell in question
%     name_list{maps,1} = mode(trace2folder(trace2folder(:,2)==map_list(maps),4));
%     %also store the cell name
%     name_list{maps,2} = uni_cells{name_list{maps,1}};    
% end
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
%% Save the maps

%define the output path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_mapmaker_out';

%assemble the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_EImaps',ori_setup,'.mat');

%save the file
save(fullfile(save_path,save_name),'all_maps')