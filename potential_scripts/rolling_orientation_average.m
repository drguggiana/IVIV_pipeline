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

% define the window width (in degrees)
window = 45;
% get the orientation and the parameter of interest
orientation_vector = cat(1,str.ORIpref);


% parameter_vector = cat(1,str.frac_vert);
% parameter_vector = sum(parameter_vector(:,3:5),2);
% parameter_vector = cat(1,str.ang_exL23);
% parameter_vector = abs(parameter_vector(:,3) - parameter_vector(:,1));
%parameter_vector = 90-abs(parameter_vector(:,5));
%parameter_vector = (abs(out_ang_exL23(:,3)-out_ang_exL23(:,1)))*69;
parameter_vector = 90-abs(out_ang_inL23(:,5));

% add osi cutoff
cutoff_osi = 0.25;
% get the osi
osi = cat(1,str.OSIpref);
% get the vector
osi_selection = osi > cutoff_osi;
% get the selected traces
orientation_vector = orientation_vector(osi_selection);
parameter_vector = parameter_vector(osi_selection);

% run the function
[rolling_orientation,rolling_ori_error] = ...
    rolling_circular_average(orientation_vector,parameter_vector,window,'ori');
%% Generate a surrogate computation for a CI

% define the number of shuffles
shuffle_number = 500;
% allocate memory to store the results
shuffle_ori = zeros(shuffle_number,180);
% take only the non-nan orientation and their values
nonnan_ori = orientation_vector(~isnan(orientation_vector));
nonnan_param = parameter_vector(~isnan(orientation_vector));

% for all the shuffles
for shuffles = 1:shuffle_number
    % randomize the parameter_vector
    rand_param = nonnan_param(randperm(length(nonnan_param)));
    % run the function
    [shuffle_ori(shuffles,:),~] = rolling_circular_average(nonnan_ori,rand_param,window,'ori');
end
% get the mean and bounds
mean_shuffle = nanmean(shuffle_ori,1);
CI_shuffle = cat(1,abs(prctile(shuffle_ori,5,1)-mean_shuffle),prctile(shuffle_ori,95,1)-mean_shuffle);
%% Plot the results

%close all
figure;set(gcf, 'Position', [800, 200, 400, 300])
set(gcf,'color','w');
shadedErrorBar(1:180,rolling_orientation,rolling_ori_error,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:180,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
xlabel('Orientation (deg)')
%ylabel('Parameter')
title(strjoin({'Rolling orientation average','window',...
    num2str(window)},'_'),'Interpreter','None')
axis tight
%% 
ylabel('L23 CoM (µm)','Color','b');
ylim([10 70]);
set(gca,'FontSize',10);
%% Calculate the rolling average for direction
% get the orientations 
% direction = cat(1,round(cat(1,str.DIRpref)),round(cat(1,str.ORIpref))+360);
direction = cat(1,str.DIRpref);
% get the exc fraction in layer 4 
% parameter = cat(1,str.frac_vert);
% parameter = sum(parameter(:,6:7),2);
parameter = cat(1,str.ang_inL23);
parameter = abs(parameter(:,3) - parameter(:,1));
% parameter_vector = 90-abs(parameter_vector(:,5));

% define the window width (in degrees)
window = 45;

% add osi cutoff
cutoff_dsi = 0.3;
% get the osi
dsi = cat(1,str.DSIpref);
% get the vector
dsi_selection = dsi > cutoff_dsi;
% get the selected traces
direction = direction(dsi_selection);
parameter = parameter(dsi_selection);

% run the function
[rolling_direction,rolling_dir_error] = ...
    rolling_circular_average(direction,parameter,window,'dir');
%% Generate a surrogate computation for a CI

% define the number of shuffles
shuffle_number = 500;
% allocate memory to store the results
shuffle_dir = zeros(shuffle_number,360);
% take only the non-nan orientation and their values
nonnan_dir = direction(~isnan(direction));
nonnan_param = parameter(~isnan(direction));

% for all the shuffles
for shuffles = 1:shuffle_number
    % randomize the parameter_vector
    rand_param = nonnan_param(randperm(length(nonnan_param)));
    % run the function
    [shuffle_dir(shuffles,:),~] = rolling_circular_average(nonnan_dir,rand_param,window,'dir');
end
% get the mean and bounds
mean_shuffle = mean(shuffle_dir,1);
CI_shuffle = cat(1,abs(prctile(shuffle_dir,5,1)-mean_shuffle),prctile(shuffle_dir,95,1)-mean_shuffle);
%% Plot the results
close all
figure

shadedErrorBar(1:360,rolling_direction,rolling_dir_error,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:360,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
xlabel('Direction')
ylabel('Parameter')
title(strjoin({'Rolling direction average','window',...
    num2str(window)},'_'),'Interpreter','None')
axis tight

