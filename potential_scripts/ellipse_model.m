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
    
    corrected_centroid = weighted_centroid(max_idx).WeightedCentroid;
    corrected_centroid(1) = (corrected_centroid(1)-8.5).*69;
    corrected_centroid(2) = (6.5-corrected_centroid(2)).*69;
    % get the centroid coordinates
    centroid_vector(maps,:) = (corrected_centroid - target_somas(maps,:))/69;
end

figure

% calculate the angle
angle_soma = rad2deg(atan2(centroid_vector(:,2),centroid_vector(:,1))+pi/2);
angle_soma(angle_soma > 180) = angle_soma(angle_soma > 180) - 360;

% angle_soma(angle_soma > 90) = angle_soma(angle_soma > 90) - 90;
% angle_soma(angle_soma < -90) = angle_soma(angle_soma < -90) + 90;

original_centroid = vertcat(str(selection_vector).ang_inL23);
original_y = abs(original_centroid(:,4) - original_centroid(:,2));
original_x = abs(original_centroid(:,3) - original_centroid(:,1));
x = original_x;
% x = abs(centroid_vector(:,1));
% x = vertcat(str(selection_vector).pialD);
% x = vertcat(str(selection_vector).ORIpref);

y = original_y;
% x = original_centroid;
% x = angle_soma;
% y = vertcat(str(selection_vector).ORIpref);
% y = abs(centroid_vector(:,2));
% c = angle_soma;
% y = vertcat(str(selection_vector).pialD);
c = vertcat(str(selection_vector).pialD);
scatter(x,y,30,c,'filled')
[rho,pval] = corr(x,y);
title(strjoin({'corr',num2str(rho),'pval',num2str(pval)},'_'),'Interpreter','None')
colormap(hsv)

% plot(vertcat(str(selection_vector).ORIpref),centroid_vector(:,1),'o')

figure
a = vertcat(str(selection_vector).ORIpref);
% a = vertcat(str(selection_vector).frac_vert);
% a = sum(a(:,3:5),2);
% a = vertcat(str(selection_vector).morph);
% a = a(:,13);
% a = vertcat(str(selection_vector).Sigmapref);
% scatter3(abs(centroid_vector(:,1)),abs(centroid_vector(:,2)),vertcat(str(selection_vector).ORIpref))
scatter3(x,y,c,[],a,'filled')
xlabel('Cx')
ylabel('Cy')
zlabel('PialD')
set(gca,'ZDir','reverse')

% define the angle range
angle_range = 0:360;
% run the calculation
[new_centroid,corr_pval] = rotate_ellipse(centroid_vector,target_ori,angle_range,angle_span,'ori');
%% Centroid test

close all
% define the number of conditions
condition_number = 5;
% define the matrix to test the different conditions
condition_matrix = zeros(5,5,condition_number);

% define manually the different conditions
% condition 1
condition_matrix([2,4],:,1) =  1;
condition_matrix(:,[2 4],1) = 0.5;
% condition 2
condition_matrix([1 2 3],[2 4],2) = 1;
condition_matrix(1,3,2) = 0.5;
% condition 3
condition_matrix([1 2 3],2,3) = 1;
condition_matrix([1 2],4,3) = 1;
condition_matrix(1,3,3) = 1;
% condition 4
condition_matrix(1,:,4) = 0.3;
condition_matrix(2,[2 3 4],4) = 0.5;
condition_matrix(3,3,4) = 1;
% condition 5
condition_matrix([1 2],[1 2],5) = 1;
condition_matrix([1 2],[4 5],5) = 0.3;
condition_matrix(1,3,5) = 0.5;

figure
% allocate memory for the centroids
centroid_matrix = zeros(condition_number,2);
% run through all the conditions and calculate the weighted centroid
for condition = 1:condition_number
    % plot the condition
    subplot(round(sqrt(condition_number)),ceil(sqrt(condition_number)),condition)
    imagesc(condition_matrix(:,:,condition))
    hold on
    % calculate the weighted centroid
    weighted_centroid = regionprops(condition_matrix(:,:,condition)~=0,...
        condition_matrix(:,:,condition),'WeightedCentroid');
    centroid_matrix(condition,:) = weighted_centroid.WeightedCentroid;
    % plot the centroid
    plot(centroid_matrix(condition,1),centroid_matrix(condition,2),'ob',...
        'MarkerFaceColor','b')
    colormap(hot)
    axis square
    
end
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