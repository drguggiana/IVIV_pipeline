function [rolling_average,rolling_std] = rolling_circular_average(circular_vector,parameter_vector,window,dir_ori_flag)
% function to take a circular rolling window average and std of the
% parameter

% define the flag dependent quantities
if strcmp(dir_ori_flag,'ori')
    single_angle = 180;
elseif strcmp(dir_ori_flag,'dir')
    single_angle = 360;
end

% get the orientations 
% orientation = cat(1,round(cat(1,str.ORIpref)),round(cat(1,str.ORIpref))+single_angle);
orientation = cat(1,round(circular_vector),round(circular_vector)+single_angle);
% get the exc fraction in layer 4 
% parameter = cat(1,str.frac_vert);
% parameter = repmat(sum(parameter_vector(:,6:7),2),2,1);
parameter = repmat(parameter_vector,2,1);

% define the window width (in degrees)
% window = 21;

% allocate memory for the average
rolling_mean_ori = zeros(2*single_angle,1);
rolling_std_ori = zeros(2*single_angle,1);
% for all degrees, starting at the window width and stopping 1 window width
% before the end
for degrees = (window-1)/2:2*single_angle-(window-1)/2
    % get the indexes of the involved cells
    indexes = orientation > degrees-(window-1)/2 &...
        orientation < degrees+(window-1)/2;
    % get the average of the involved numbers
    rolling_mean_ori(degrees) = nanmean(parameter(indexes));
    rolling_std_ori(degrees) = nanstd(parameter(indexes))./sqrt(sum(indexes));
end

% assemble the orientation only vector based on the window width
rolling_average = cat(1,rolling_mean_ori(single_angle+1:single_angle-1+(window-1)/2),...
    rolling_mean_ori((window-1)/2:single_angle));
rolling_std = cat(1,rolling_std_ori(single_angle+1:single_angle-1+(window-1)/2),...
    rolling_std_ori((window-1)/2:single_angle));