function bsub_map = trace_fetch(map_path,trace_range,trace_background,num_positions)
      
%load the xsg file
temp_xsg = load(map_path,'-mat');
%check the number of grid points in the cell
temp_map = temp_xsg.header.mapper.mapper.mapPatternArray;

%load the time trace
temp_trace = temp_xsg.data.ephys.trace_1;
%reshape to trim the relevant part of the data only
temp_trace = reshape(temp_trace,[],num_positions);

%get the background activity
background_act = mean(temp_trace(trace_background,:),1);

%trim the trace
temp_trace = temp_trace(trace_range,:);
%subtract background activity
temp_trace = temp_trace - background_act;
%also load the map order to have all the maps ordered in their location
%instead of in time of stimulation
%linearize the map and apply to the data
bsub_map = temp_trace(:,temp_map(:));