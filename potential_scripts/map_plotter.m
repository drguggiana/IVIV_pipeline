function map_matrix = map_plotter(map_info,map_data)
% plot the input map based on the input traces

% %save the cell's name
% all_maps{cells,1} = uni_cells{cells};
% %also the cells soma coordinates
% all_maps{cells,8} = soma_unique(cells,:);

% %get the info for the traces from this cell
% curr_cell = trace2folder(trace2folder(:,4)==cells,:);
%

% allocate memory to output the maps
map_matrix = cell(3,1);
%for the polarities present
for polars_count = 1:3
%     figure
    
    
    if polars_count > 2
        polars = 2;
    else
        polars = polars_count;
    end
    %
    % get this map's info
    curr_traces = map_info{polars};
    % get the current map
    bsub_map = map_data(:,:,polars);
    % get the number of pixels
    num_pixels = size(bsub_map,2);
    % allocate memory for the final map
    proc_maps = zeros(sqrt(num_pixels),sqrt(num_pixels),1);
    %
    %     %get the set of maps present in this cell and polarity
    %     curr_maps = unique(curr_cell(curr_cell(:,3)==polars-1,2));
    %     %get the paths to the cells maps
    %     cell_paths = folder_all(curr_maps);
    %     %allocate memory to store the processed maps
    %     proc_maps = zeros(num_positions,length(cell_paths));
    %     %allocate memory to store the 16x16x6 trace2folder array
    %     pixel_sub = zeros(num_positions,size(trace2folder,2),length(cell_paths));
    %
    %
    %
    %     %for all the maps
    %     for maps = 1:length(cell_paths)
    %         %also get the info for the current traces
    %         curr_traces = curr_cell(curr_cell(:,3)==polars-1&curr_cell(:,2)==curr_maps(maps),:);
    %         %if there are interpolated traces, load the corresponding
    %         %map
    %         if exist('interp_cell','var')
    %             %define the target (synaptic) window
    %             target_window = 8:157;
    %             %load the map from the interpolation cell
    %             bsub_map = interp_cell{curr_maps(maps)};
    %         else
    %             %define the target (synaptic) window
%     target_window = 71:1570;
    target_window = 1:150;
    %             %fetch the background subtracted map from the raw data
    %             bsub_map = trace_fetch(cell_paths{maps},trace_range,trace_background,num_positions);
    %         end
    %         %count the types of responses for this map
    %         map_info = response_counter(curr_traces,polars,bsub_map);
    %         %store the trace2folder info
    %         pixel_sub(:,:,maps) = curr_traces;
    
    %                     %if there is a single rep for the map, skip it
    %                     if length(cell_paths) < 2
    %
    %                          proc_maps(:,maps) = zeros(num_positions,1);
    %
    %                     else
    
    %process the traces according to the type of response
    %for all the traces
    for trace = 1:num_pixels
        %get the response type
        resp_type = curr_traces(trace,5);
        %if it's a NaN, skip the processing (i.e. blank trace with a zero)
        if isnan(resp_type)
            continue
        end
        % get the onset
        onset = curr_traces(trace,7);
        % get the fs flag based on the trace onset
        if polars_count == 2
            if onset > 100
                fs_flag = 1;
            else
                fs_flag = 0;
            end
        else
            if onset < 100
                fs_flag = 1;
            else
                fs_flag = 0;
            end
        end
        %process the trace accordingly
        proc_maps(trace) = Trace_process(bsub_map(target_window,trace),resp_type,polars_count,fs_flag);
        
    end
    
%     imagesc(proc_maps);
    
    map_matrix{polars_count} = proc_maps;
    %                     end
    
end

% concatenate the maps
map_matrix = cat(3,map_matrix{:});

%
%
%     %blank the positions that only have a map in a single repetition
%     %get the positions
%     single_blank = prod(proc_maps,2)==0;
%     proc_maps(single_blank,:) = 0;
%
% %     %save the average maps
% %     all_maps{cells,polars+1} = nanmean(reshape(proc_maps,sqrt(num_positions),sqrt(num_positions),[]),3);
% %     %store the extra info
% %     all_maps{cells,polars+3} = map_info;
% %     %and also the pixel by pixel info in an array
% %     all_maps{cells,polars+5} = reshape(mode(pixel_sub,3),sqrt(num_positions),sqrt(num_positions),[]);
% end