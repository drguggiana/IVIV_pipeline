%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Calculate the rolling average for orientation

% get the orientations 
orientation = cat(1,round(cat(1,str.ORIpref)),round(cat(1,str.ORIpref))+180);
% get the exc fraction in layer 4 
parameter = cat(1,str.frac_vert);
parameter = repmat(sum(parameter(:,23:24),2),2,1);

% define the window width (in degrees)
window = 21;

% allocate memory for the average
rolling_mean_ori = zeros(360,1);
rolling_std_ori = zeros(360,1);
% for all degrees, starting at the window width and stopping 1 window width
% before the end
for degrees = (window-1)/2:360-(window-1)/2
    % get the indexes of the involved cells
    indexes = orientation > degrees-(window-1)/2 &...
        orientation < degrees+(window-1)/2;
    % get the average of the involved numbers
    rolling_mean_ori(degrees) = mean(parameter(indexes));
    rolling_std_ori(degrees) = std(parameter(indexes))./sqrt(sum(indexes));
end

% assemble the orientation only vector based on the window width
rolling_orientation = cat(1,rolling_mean_ori(181:180-1+(window-1)/2),...
    rolling_mean_ori((window-1)/2:180));
rolling_ori_error = cat(1,rolling_std_ori(181:180-1+(window-1)/2),...
    rolling_std_ori((window-1)/2:180));
%% Plot the results
close all
figure

shadedErrorBar(1:180,rolling_orientation,rolling_ori_error)
xlabel('Orientation')
ylabel('Parameter')
title(strjoin({'Rolling orientation average','window',...
    num2str(window)},'_'),'Interpreter','None')
%% Calculate the rolling average for direction
% get the orientations 
direction = cat(1,round(cat(1,str.DIRpref)),round(cat(1,str.ORIpref))+360);
% get the exc fraction in layer 4 
parameter = cat(1,str.frac_vert);
parameter = repmat(sum(parameter(:,23:24),2),2,1);

% define the window width (in degrees)
window = 21;

% allocate memory for the average
rolling_mean_dir = zeros(720,1);
rolling_std_dir = zeros(720,1);
% for all degrees, starting at the window width and stopping 1 window width
% before the end
for degrees = (window-1)/2:720-(window-1)/2
    % get the indexes of the involved cells
    indexes = orientation > degrees-(window-1)/2 &...
        orientation < degrees+(window-1)/2;
    % get the average of the involved numbers
    rolling_mean_dir(degrees) = mean(parameter(indexes));
    rolling_std_dir(degrees) = std(parameter(indexes))./sqrt(sum(indexes));
end

% assemble the orientation only vector based on the window width
rolling_direction = cat(1,rolling_mean_dir(361:360-1+(window-1)/2),...
    rolling_mean_dir((window-1)/2:360));
rolling_dir_error = cat(1,rolling_std_dir(361:360-1+(window-1)/2),...
    rolling_std_dir((window-1)/2:360));
%% Plot the results
close all
figure

shadedErrorBar(1:360,rolling_direction,rolling_dir_error)
xlabel('Orientation')
ylabel('Parameter')
title(strjoin({'Rolling orientation average','window',...
    num2str(window)},'_'),'Interpreter','None')
