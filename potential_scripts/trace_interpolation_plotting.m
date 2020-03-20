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