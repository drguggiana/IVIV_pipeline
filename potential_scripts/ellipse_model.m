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
% y = vertcat(str(selection_vector).ORIpref);

y = original_y;
% x = original_centroid;
% x = angle_soma;
% y = vertcat(str(selection_vector).ORIpref);
% y = abs(centroid_vector(:,2));
% c = angle_soma;
% y = vertcat(str(selection_vector).pialD);
c = vertcat(str(selection_vector).ORIpref);

frac = vertcat(str(selection_vector).frac_vert);
frac = sum(frac(:,6:7),2);

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
%% Linear regression of Cx and Cy
close all

% X = abs(reshape(target_maps,size(target_maps,1)*size(target_maps,2),number_of_targets)');
% X = [x,y(randperm(length(y)))];
% X = sqrt(x.^2+y.^2);
X = [x,y];
X = normr_2(X,2);
X(isnan(X)) = 0;

y_m = vertcat(str(selection_vector).ORIpref);
% y = normr_2(centroid_vector(:,1));
mdl = fitlm(X,y_m);

figure
plot(mdl)

mdl_svm = fitrsvm(X,y_m,'KernelFunction','gaussian');
y_pred = predict(mdl_svm,X);
figure
scatter(X,y_pred)
hold on
scatter(X,y_m)
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
%% Attempt with cutting the cloud
close all

% define a vector of distance to try
distance_vector = 10:10:200;
distance_vector = horzcat(distance_vector,8*69);

distance_number = length(distance_vector);
% allocate memory to store the results
distance_correlation = zeros(distance_number,2);
% for all the distances
for distance = 1:distance_number
    
    fprintf(strjoin({'Current distance:',num2str(distance),...
        'out of',num2str(distance_number),'\r\n'},'_'))
    % define the distance at which to cut
    distance_threshold = distance_vector(distance);
    
    % define the target angle for selecting the maps
    target_angle = 125;
    % define angle_span
    angle_span = 25;
    
    % get the maps
%     selection_vector = vertcat(str.OSIpref)>0.25;
%     selection_vector = selection_vector & (vertcat(str.ORIpref)<(angle_span+target_angle)...
%         & vertcat(str.ORIpref)>(angle_span-target_angle));
    selection_vector = vertcat(str.resp)==0;
    
    target_maps = cat(3,str(selection_vector).subpixel_inhMap);
    % get the number of maps
    target_num = size(target_maps,3);
    % also the soma positions
    target_somas = cat(1,str(selection_vector).subpixel_soma);
    % also the orientations
    target_ori = vertcat(str(selection_vector).ORIpref);
    % get the original centroids
    original_centroids = vertcat(str(selection_vector).ang_inL23);
    cx = abs(original_centroids(:,3)-original_centroids(:,1));
    cy = abs(original_centroids(:,4)-original_centroids(:,2));
    % calculate the length
    original_centroids = sqrt(cx.^2 + cy.^2);
    
    % get the number of maps
    number_of_targets = size(target_maps,3);
    
    % % select a layer
    % target_maps = target_maps(3:5,:,:);
    
    % define the grid spacing in microns
    grid_spacing = 69;
    
    % get the map size in microns
    map_size = round(16.*grid_spacing);
    % get the map limits
    map_lim = map_size/2-grid_spacing/2;
    
    % create the grid
    [Y,X] = ndgrid(-map_lim:grid_spacing:map_lim,...
        -map_lim:grid_spacing:map_lim);
    
    % define the single micron grid
    [Y_single,X_single] = ndgrid(-map_lim:map_lim,...
        -map_lim:map_lim);
    
    % get the slice ori cloud centers for setup 2
    centroid_vector = zeros(target_num,1);
    
    % for all the cells
    for cells = 1:target_num
        %     % get the soma centers in x
        %     soma = str(cells).somaCenter(1);
        %     % flip the sign if it's a sliceOri 0 map
        %     if str(cells).sliceOri == 0
        %         soma = -soma;
        %     end
        %     % generate the centered grid
        %     center_X = X - soma;
        
        % get the corresponding map
        map = target_maps(:,:,cells);
        % interpolate the map
        interpolant = griddedInterpolant(Y,X,map);
        interp_map = interpolant(Y_single,X_single);
        % determine the cutting position
        cutting_location = target_somas(cells,1) + distance_threshold;
        % find the position in the actual map
        [~,cutting_map] = min(abs(X_single(1,:)-cutting_location));
        % cut the cloud
        cut_map = interp_map;
        cut_map(:,cutting_map:end) = 0;
        
%         cutting_location2 = target_somas(cells,1) - distance_threshold;
%         % find the position in the actual map
%         [~,cutting_map2] = min(abs(X_single(1,:)-cutting_location2));
% %         % cut the cloud
% %         cut_map = interp_map;
%         % cut the map symmetrically
%         cut_map(1:cutting_map:end) = 0;
        % select the layer for centroid calculation
        cut_map = cut_map(3*grid_spacing:5*grid_spacing,:);
        % get the components
        cc = bwconncomp(cut_map);
        rp = regionprops(cc,cut_map,{'Area','WeightedCentroid'});
        % if there's more than one, leave the largest
        if cc.NumObjects > 1
            areas = cat(1,rp.Area);
            [~,idx] = max(areas);
            rp = rp(idx);
        end
        % store the centroid
        centroid_coord = rp.WeightedCentroid;
        cx = (X_single(1,round(centroid_coord(1))) - target_somas(cells,1))./69;
        cy = (-Y_single(round(centroid_coord(2)),1) - 2*69 - target_somas(cells,2))./69;
        % calculate the length and store
        centroid_vector(cells) = sqrt(cx^2+cy^2);
    end
    % calculate and store the correlation and pvalue
    [distance_correlation(distance,1),distance_correlation(distance,2)] =...
        corr(centroid_vector,original_centroids);
end
% % plot the 2
% figure
% scatter(original_centroids,target_ori)
% hold on
% scatter(centroid_vector,target_ori)
% % scatter(abs(centroid_vector(:,1)-target_somas(:,1))./69,target_ori)
% ylabel('Orientation')
% xlabel('Asymmetry')
% [rho1,pval1] = corr(original_centroids,target_ori);
% [rho2,pval2] = corr(centroid_vector,target_ori);
% title(strjoin({num2str(rho1),num2str(pval1),num2str(rho2),num2str(pval2)},'_'),'Interpreter','None')
%% plot the results

close all

figure
plot(distance_vector,distance_correlation(:,1))
hold on
plot(distance_vector,distance_correlation(:,2))

plot([distance_vector(1),distance_vector(end)],[0.05 0.05],'r')
axis tight
xlabel('Cutting distance from soma')
[rho1,pval1] = corr(original_centroids,target_ori);

title(strjoin({num2str(rho1),num2str(pval1)},'_'),'Interpreter','None')
legend({'Correlation','p value'},'location','northwest')
