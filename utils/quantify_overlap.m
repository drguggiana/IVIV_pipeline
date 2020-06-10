function [overlap_index,pairs_cell,overlap_matrix] = quantify_overlap(idx1,idx2,noise_threshold)
% quantify the maximum possible overlap between 2 sets of clusters

% get the list and number of clusters for each set
clu_list1 = unique(idx1);
clu_num1 = length(clu_list1);

clu_list2 = unique(idx2);
clu_num2 = length(clu_list2);

% calculate the actual noise threshold
noise_number = ceil(noise_threshold*length(idx1));

% get the list and number of possible pairwise combinations
% combo_list = nchoosek(1:(clu_num1+clu_num2),2);
[combos_1,combos_2] = ndgrid(1:clu_num1,1:clu_num2);
% flatten the lists
combos_1 = combos_1(:);
combos_2 = combos_2(:);
combo_num = size(combos_1,1);

% allocate memory for the results
overlap_pairs = zeros(clu_num1,clu_num2);

% for all the combinations
for combo = 1:combo_num
    % calculate the overlap
    overlap_pairs(combos_1(combo),combos_2(combo)) = ...
        sum((idx1==combos_1(combo))&(idx2==combos_2(combo)));
end

% copy the overlap matrix and output
overlap_matrix = overlap_pairs;
% get the minimum number of clusters
min_clu = min(clu_num1,clu_num2);

% allocate memory for the max scores
scores_vector = zeros(min_clu,1);
% allocate memory for the pairs
pairs_cell = cell(max(clu_num1,clu_num2),2);

% for the min number of clusters
for clu = 1:min_clu
    % get the marginal sums
    marginal_x = sum(overlap_pairs,2);
    marginal_y = sum(overlap_pairs,1);
    % once all numbers are below the noise_number, reset to 0
    if max(overlap_pairs(:))<=noise_number
        noise_number = 0;
    end
    
    % detect which one has the highest number
    if max(marginal_x) > max(marginal_y)
        % find row with the max
        [max_val,max_row] = max(marginal_x);
        % identify the rows involved
        idx_cols = overlap_pairs(max_row,:)>noise_number;
        % save the max and the pairings
        scores_vector(max_row) = max_val;
        pairs_cell{max_row,1} = find(idx_cols);
        % blank the involved rows and columns
        overlap_pairs(max_row,:) = 0;
        overlap_pairs(:,idx_cols) = 0;
    else
        % find row with the max
        [max_val,max_col] = max(marginal_y);
        % identify the rows involved
        idx_rows = overlap_pairs(:,max_col)>noise_number;
        % save the max and the pairings
        scores_vector(max_col) = max_val;
        pairs_cell{max_col,2} = find(idx_rows);
        % blank the involved rows and columns
        overlap_pairs(:,max_col) = 0;
        overlap_pairs(idx_rows,:) = 0;
    end
end

% calculate the overlap index
overlap_index = sum(scores_vector)/length(idx1);