clearvars
close all
Paths
% define the path to the structure
% str_path = 'R:\Share\Simon\Drago_Volker_Simon\Full_data_structure\str_iviv_200323.mat';
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
%% Regress a linear model on the data
X = abs(reshape(target_maps,size(target_maps,1)*size(target_maps,2),number_of_targets)');
X = normr_2(X,2);
X(isnan(X)) = 0;
y = normr_2(centroid_vector(:,1));
mdl = fitlm(X,y);
%% Plot the results
close all

% get the coefficients
coefficients = mdl.Coefficients.Estimate;
% exclude the intercept and reshape to map
coefficients = reshape(coefficients(2:end),size(target_maps,1),size(target_maps,2));
% plot the map
figure
imagesc(log(abs(coefficients)))
%% Use SVR for the fit

% X = abs(reshape(target_maps,size(target_maps,1)*size(target_maps,2),number_of_targets)');
% X = normr_2(X,2);
% X(isnan(X)) = 0;
% y = normr_2(centroid_vector(:,1));
mdl_svm = fitrsvm(X,y);
%% Plot the SVM results

close all

% get the coefficients
coefficients = mdl_svm.Beta;
% exclude the intercept and reshape to map
coefficients = reshape(coefficients,size(target_maps,1),size(target_maps,2));
% plot the map
figure
imagesc((coefficients))