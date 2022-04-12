function [conf_pearson conf_spearman conf_kendall] = conf_interval(variable1, variable2)
% define the target variables (observations in rows, features in columns)
% variable1 = morpho_props;
% variable2 = pial_depth;
% define the alpha
alpha = 0.05;
% define the number of shuffles
num_shuffle = 1000;
% allocate memory for the output
shuffle_mat = zeros(3,num_shuffle,size(variable1,2),size(variable2,2));
% exclude rows with nans in the variables
nonan_vector = ~any(isnan(variable1),2) & ~any(isnan(variable2),2);
variable1 = variable1(nonan_vector, :);
variable2 = variable2(nonan_vector, :);
% get the number of features in var1
num_variable1 = size(variable1,2);
% for all the shuffles
for shuffle = 1:num_shuffle
    % get the random sample with replacement
    random_sample = datasample([variable1,variable2],size(variable1,1));
    shuffle1 = random_sample(:,1:num_variable1);
    shuffle2 = random_sample(:,num_variable1+1:end);
    [rho_pearson_shuffle,~] = corr(shuffle1, shuffle2, 'type', 'Pearson', 'rows', 'pairwise');
    [rho_spearman_shuffle,~] = corr(shuffle1, shuffle2, 'type', 'Spearman', 'rows', 'pairwise');
    [rho_kendall_shuffle,~] = corr(shuffle1, shuffle2, 'type', 'Kendall', 'rows', 'pairwise');
    % store the correlation matrices
    shuffle_mat(1,shuffle,:,:) = rho_pearson_shuffle;
    shuffle_mat(2,shuffle,:,:) = rho_spearman_shuffle;
    shuffle_mat(3,shuffle,:,:) = rho_kendall_shuffle;
end
% calculate the confidence limits
confidence_up = 100-(alpha*100/2);
confidence_down = alpha*100/2;
% calculate the confidence intervals for each one
conf_pearson = [prctile(squeeze(shuffle_mat(1,:,:,:)),confidence_down,1);...
    prctile(squeeze(shuffle_mat(1,:,:,:)),confidence_up,1)];
conf_spearman = [prctile(squeeze(shuffle_mat(2,:,:,:)),confidence_down,1);...
    prctile(squeeze(shuffle_mat(2,:,:,:)),confidence_up,1)];
conf_kendall = [prctile(squeeze(shuffle_mat(3,:,:,:)),confidence_down,1);...
    prctile(squeeze(shuffle_mat(3,:,:,:)),confidence_up,1)];
end