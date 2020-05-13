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