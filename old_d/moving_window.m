function data_out = moving_window(data_in,window_size,function_handle,position_flag)
%this function applies the function "function_handle" in a rolling window
%across the first dimension of the data, with window size determined by
%window_size
%position_flag controls whether the window is centered on the point (0) or it
%precedes it (1)

%get the size of the data
data_num = size(data_in,1);

%allocate memory for the output
data_out = zeros(size(data_in));

%for the moving window range
for window = 1:data_num
    
    switch position_flag
        case 0
            %define the behavior depending on the position
            if window < ceil(window_size/2)
                %get the actual window
                window_idx = 1:window_size;
            elseif window > data_num - floor(window_size/2)
                %get the actual window
                window_idx = data_num-window_size:data_num;
            else
                %get the actual window
                window_idx = window-round(window_size/2)+1:window+round(window_size/2)-1;
            end
        case 1
            %define the behavior depending on the position
            if window < window_size
                %get the actual window
                window_idx = 1:window_size;
            else
                %get the actual window
                window_idx = window-window_size+1:window;
            end
    end
    
    %calculate the value
    data_out(window,:) = function_handle(data_in(window_idx,:));
    
end

% figure
% plot(data_in(:,1))
% hold on
% plot(data_out(:,1))