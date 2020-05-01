clearvars
close all
Paths
% define the path to the structure
str_path = structure_file_path;
% load the structure
load(str_path)
%% Calculate the rotated centroids
close all


% fit an ellipse to the input clouds from the maps that are aligned along
% the slice

% get the maps 
selection_vector = vertcat(str.OSIpref)>0.25;
% selection_vector = ones(cell_num,1)==1;
target_maps = cat(3,str(selection_vector).subpixel_inhMap);
% also the soma positions
target_somas = cat(1,str(selection_vector).subpixel_soma);
% also the orientations
target_ori = vertcat(str(selection_vector).ORIpref);

% get the number of maps
number_of_targets = size(target_maps,3);

% select a layer
target_maps = target_maps(3:5,:,:);
% define angle_span
angle_span = 25;

% plot it
% figure
% imagesc(target_average)

% allocate memory to store the centroids
centroid_vector = zeros(number_of_targets,2);
% for all the maps
for maps = 1:number_of_targets
    % calculate the position of the weighted centroid
    weighted_centroid = regionprops(target_maps(:,:,maps)~=0,target_maps(:,:,maps),'Area','WeightedCentroid');
    if isempty(weighted_centroid)
        centroid_vector(maps,:) = [NaN,NaN];
        continue
    end
    %     weighted_centroid = regionprops(target_maps(:,:,maps)~=0,'Area','Centroid');
    % get the largest CC by area
    target_areas = cat(1,weighted_centroid.Area);
    [~,max_idx] = max(target_areas);
    % get the centroid coordinates
    centroid_vector(maps,:) = abs((weighted_centroid(max_idx).WeightedCentroid-8.5).*69 - target_somas(maps,1));
end

% define the angle range
angle_range = 0:360;
% run the calculation
[new_centroid,corr_pval] = rotate_ellipse(centroid_vector,target_ori,angle_range,angle_span,'ori');
%% Shuffle the angles to obtain a CI

% define the number of shuffles
shuffle_number = 100;

% allocate memory for the results
shuffle_results = zeros(shuffle_number,length(angle_range),2);

% for all of the shuffles
for shuffle = 1:shuffle_number
    % shuffle the angles
    shuffle_ori = target_ori(randperm(number_of_targets));
    % run the calculation
    [~,temp_corr] = rotate_ellipse(centroid_vector,shuffle_ori,angle_range,angle_span,'ori');
    % store the results
    shuffle_results(shuffle,:,:) = temp_corr;
end

% get the confidence intervals
mean_shuffle = nanmean(shuffle_results,1); 
shuffle_CI = cat(1,abs(prctile(shuffle_results,5,1)-mean_shuffle),...
    prctile(shuffle_results,95,1)-mean_shuffle);
%% Plot the results


% figure
% scatter(centroid_vector(:,1),centroid_vector(:,2),10,'ko')
% hold on
% scatter(new_centroid(1:50),centroid_vector(:,2),10,'ro')

nanmean(centroid_vector(:,1))
nanstd(centroid_vector(:,1))

figure
plot(angle_range,corr_pval(:,1),'r')
hold on
plot(angle_range,corr_pval(:,2),'b')
shadedErrorBar(angle_range,mean_shuffle(1,:,1),shuffle_CI(:,:,1),'transparent',1,'lineprops','k')
legend({'correlation','p value','shuffle +/- 95% CI'})
xlabel('Angle of alignment')