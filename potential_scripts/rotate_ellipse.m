function [new_centroid,corr_pval] = rotate_ellipse(centroid_vector,target_ori,angle_range,angle_span,oridir_flag)

switch oridir_flag
    case 'ori'
        oridir_angle = 180;
    case 'dir'
        oridir_angle = 360;
end

% duplicate centroid_vector and target_ori
centroid_vector = repmat(centroid_vector,2,1);
target_ori = [target_ori;target_ori+oridir_angle];
% allocate memory to store the results
corr_pval = zeros(length(angle_range),2);
% initialize a counter
counter = 1;

% for all possible target angles
for target_angle = angle_range
    
    % get the average angle around the current angle center
    long_average = nanmean(centroid_vector(...
        target_ori<(target_angle+angle_span)&target_ori>(target_angle-angle_span),1));
    
    % assume that as the average length of all maps and calculate the resulting
    % centroid for all maps
    % allocate memory for the new centroids
    new_centroid = zeros(size(target_ori,1),1);
    % for all the maps
    for maps = 1:size(target_ori,1)
        new_centroid(maps) = abs(cos(deg2rad(target_angle-target_ori(maps)))*long_average);
    end
    
    % eliminate nan values
    not_nan = ~isnan(centroid_vector(:,1))&~isnan(new_centroid);
    new_centroid = new_centroid(not_nan);
    % if it's empty, skip
    if isempty(new_centroid)
        % update the counter
        counter = counter + 1;
        continue
    end
    centroid_vector_nan = centroid_vector(not_nan,:);
    % run the correlation
    [corr_pval(counter,1),corr_pval(counter,2)]= corr(centroid_vector_nan(:,1),new_centroid);
    % update the counter
    counter = counter + 1;
end