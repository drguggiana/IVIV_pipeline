function rolling_avg_display(str,parameter_vector,parameter_vector2,window,or,multi)
%% Calculate the rolling average for orientation
if or==1
% define the window width (in degrees)
%window = 45;
% get the orientation and the parameter of interest
orientation_vector = cat(1,str.ORIpref);

% temp=[str(:).ORIpref]+90
% temp(temp>179) = temp(temp>179) -180;
% orientation_vector=temp';
% add osi cutoff
cutoff_osi = 0.25;
% get the osi
osi = cat(1,str.OSIpref);
% get the vector
osi_selection = osi > cutoff_osi;
% get the selected traces
orientation_vector = orientation_vector(osi_selection);
parameter_vector = parameter_vector(osi_selection);
parameter_vector2=parameter_vector2(osi_selection);

% run the function
[rolling_orientation,rolling_ori_error] = ...
    rolling_circular_average(orientation_vector,parameter_vector,window,'ori');
% run the function
[rolling_orientation2,rolling_ori_error2] = ...
    rolling_circular_average(orientation_vector,parameter_vector2,window,'ori');
%% Generate a surrogate computation for a CI

% define the number of shuffles
shuffle_number = 500;
% allocate memory to store the results
shuffle_ori = zeros(shuffle_number,180);
% take only the non-nan orientation and their values
nonnan_ori = orientation_vector(~isnan(orientation_vector));
nonnan_param = parameter_vector(~isnan(orientation_vector));
nonnan_param2 = parameter_vector2(~isnan(orientation_vector));

% for all the shuffles
for shuffles = 1:shuffle_number
    % randomize the parameter_vector
    rand_param = nonnan_param(randperm(length(nonnan_param)));
    rand_param2 = nonnan_param2(randperm(length(nonnan_param2)));
    % run the function
    [shuffle_ori(shuffles,:),~] = rolling_circular_average(nonnan_ori,rand_param,window,'ori');
    [shuffle_ori2(shuffles,:),~] = rolling_circular_average(nonnan_ori,rand_param2,window,'ori');
end
% get the mean and bounds
mean_shuffle = nanmean(shuffle_ori,1);
mean_shuffle2 = nanmean(shuffle_ori2,1);
CI_shuffle = cat(1,abs(prctile(shuffle_ori,5,1)-mean_shuffle),prctile(shuffle_ori,95,1)-mean_shuffle);
CI_shuffle2 = cat(1,abs(prctile(shuffle_ori2,5,1)-mean_shuffle2),prctile(shuffle_ori2,95,1)-mean_shuffle2);
%% Plot the results
if multi==1
%close all
figure;set(gcf, 'Position', [800, 500, 200, 225])
set(gcf,'color','w');
shadedErrorBar(1:180,rolling_orientation,rolling_ori_error,'transparent',1,'lineprops','r')
hold on
shadedErrorBar(1:180,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
hold on
shadedErrorBar(1:180,rolling_orientation2,rolling_ori_error2,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:180,mean_shuffle2,CI_shuffle2,'transparent',1,'lineprops','k')
xlabel('Orientation (deg)');
xlim([0 180]);xticks([0:45:180]);
lgd=legend('data','shuffle');legend boxoff;lgd.Location='northwest';
lgd.ItemTokenSize = [5,5];
%ylabel('Parameter')
%title(strjoin({'Rolling orientation average','window',...
   % num2str(window)},'_'),'Interpreter','None')
axis tight
else multi==0
    figure;set(gcf, 'Position', [800, 500, 200, 225])
set(gcf,'color','w');
shadedErrorBar(1:180,rolling_orientation,rolling_ori_error,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:180,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
hold on
xlabel('Orientation (deg)');
xlim([0 180]);xticks([0:45:180]);
lgd=legend('data','shuffle');legend boxoff;lgd.Location='northwest';
lgd.ItemTokenSize = [5,5];
%ylabel('Parameter')
%title(strjoin({'Rolling orientation average','window',...
   % num2str(window)},'_'),'Interpreter','None')
axis tight
end
else or==2
    %window = 65;
direction_vector = cat(1,str.DIRpref);
% add osi cutoff
cutoff_dsi = 0.25;
% get the osi
dsi = cat(1,str.DSIpref);
% get the vector
dsi_selection = dsi > cutoff_dsi;
% get the selected traces
direction_vector = direction_vector(dsi_selection);
parameter_vector = parameter_vector(dsi_selection);
parameter_vector2=parameter_vector2(dsi_selection);
%orientation_vector = orientation_vector(osi_selection);
%parameter_vector = parameter_vector(osi_selection);
% run the function
[rolling_direction,rolling_dir_error] = ...
    rolling_circular_average(direction_vector,parameter_vector,window,'dir');
% run the function
[rolling_direction2,rolling_dir_error2] = ...
    rolling_circular_average(direction_vector,parameter_vector2,window,'dir');
%% Generate a surrogate computation for a CI

% define the number of shuffles
shuffle_number = 500;
% allocate memory to store the results
shuffle_dir = zeros(shuffle_number,360);
% take only the non-nan orientation and their values
nonnan_dir = direction_vector(~isnan(direction_vector));
nonnan_param = parameter_vector(~isnan(direction_vector));
nonnan_param2 = parameter_vector2(~isnan(direction_vector));
% for all the shuffles
for shuffles = 1:shuffle_number
    % randomize the parameter_vector
    rand_param = nonnan_param(randperm(length(nonnan_param)));
    rand_param2 = nonnan_param2(randperm(length(nonnan_param2)));
    % run the function
    [shuffle_dir(shuffles,:),~] = rolling_circular_average(nonnan_dir,rand_param,window,'dir');
    [shuffle_dir2(shuffles,:),~] = rolling_circular_average(nonnan_dir,rand_param2,window,'dir');
end
% get the mean and bounds
mean_shuffle = mean(shuffle_dir,1);
mean_shuffle2 = mean(shuffle_dir2,1);
CI_shuffle = cat(1,abs(prctile(shuffle_dir,5,1)-mean_shuffle),prctile(shuffle_dir,95,1)-mean_shuffle);
CI_shuffle2 = cat(1,abs(prctile(shuffle_dir2,5,1)-mean_shuffle2),prctile(shuffle_dir2,95,1)-mean_shuffle2);
%% Plot the results
if multi==1

figure;set(gcf, 'Position', [800, 500, 200, 225])
set(gcf,'color','w');
shadedErrorBar(1:360,rolling_direction,rolling_dir_error,'transparent',1,'lineprops','r')
hold on
shadedErrorBar(1:360,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
hold on;
shadedErrorBar(1:360,rolling_direction2,rolling_dir_error2,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:360,mean_shuffle2,CI_shuffle2,'transparent',1,'lineprops','k')
xlabel('Direction')
ylabel('Parameter')
title(strjoin({'Rolling direction average','window',...
    num2str(window)},'_'),'Interpreter','None')
axis tight
else
  figure;set(gcf, 'Position', [800, 500, 200, 225])
set(gcf,'color','w');
shadedErrorBar(1:360,rolling_direction,rolling_dir_error,'transparent',1,'lineprops','b')
hold on
shadedErrorBar(1:360,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')  
xlabel('Direction')
ylabel('Parameter')
title(strjoin({'Rolling direction average','window',...
    num2str(window)},'_'),'Interpreter','None')
axis tight
end
end
end