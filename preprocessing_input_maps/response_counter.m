function data_out = response_counter(curr_traces,polars,bsub_map)

%allocate memory for the output
data_out = zeros(5,1);
switch polars
    case 1 %excitation
        %get the minimum deflection in the map
        data_out(1) = min(bsub_map(:));
        %calculate the number of direct traces
        data_out(2) = sum(curr_traces(:,5)==0);
        %calculate the number of second window traces
        data_out(3) = sum(curr_traces(:,5)==2);
        %calculate the number of synaptic traces
        data_out(4) = sum(curr_traces(:,5)==1);
        
    case 2 %inhibition
        %get only the analyzed traces
        temp_trace = bsub_map(:,curr_traces(:,5)>0);
        %calculate the percentage of positive charge in the map wrt the
        %total
        data_out(1) = sum(temp_trace(temp_trace(:)>0))/sum(abs(temp_trace(:))); 
        %calculate the number of direct traces and neg bigger
        data_out(2) = sum(curr_traces(:,5)==0);
        %calculate the number of synaptic traces and pos bigger
        data_out(3) = sum(curr_traces(:,5)==1);
        %calculate the number of synaptic traces and neg bigger
        data_out(4) = sum(curr_traces(:,5)==2);
        %calculate the number of direct traces and pos bigger
        data_out(5) = sum(curr_traces(:,5)==3);
end

