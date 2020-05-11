% Model the contribution of different variables in explaining the
% orientation preference

%% Clean up and load the data

clearvars
close all
Paths
% define the path to the structure
str_path = structure_file_path;
% load the structure
load(str_path)
%% Calculate the fit with different parameters

% get the selection vector (based on the OS cells)
selection_vector = vertcat(str.OSIpref)>0.25;
% define the target variables
response = 'ORIpref';

targets = {'ang_inL23_Vt','ang_inL4_Cx','pialD','frac_vert_exL4','frac_vert_exL23',...
    'ang_exL23_Vt','ang_inL4_Vt','frac_vert_inL4','frac_vert_inL23'};

% get the number of parameters
num_parameters = length(targets);
% get the combinations of parameters
comb = nchoosek(1:num_parameters,2);
% add the singles
comb = vertcat(horzcat((1:num_parameters)',zeros(num_parameters,1)),comb);
% get the number of combinations
num_comb = size(comb,1);

% allocate memory for the saving
combination_results = zeros(num_comb,1);
% allocate a cell for the coefficients
combination_coefficients = cell(num_comb,1);
% allocate memory for the rmse
combination_rmse = zeros(num_comb,1);
% allocate memory for the labels
label_cell = cell(num_comb,1);
% for all the combinations
for combos = 1:num_comb
    % define the sub targets based on the combo
    % get the indexing vector
    idx_vector = comb(combos,:);
    idx_vector = idx_vector(idx_vector>0);
    sub_targets = targets(idx_vector);
    % run the regression
    [combination_results(combos),combination_coefficients{combos},~,~,combination_rmse(combos)] =...
        GLM_fitting(str,selection_vector,response,...
        sub_targets,0,true);
    % write the label also
    label_cell{combos} = strjoin(sub_targets,',');
end
%% Plot the parameter results

close all

% plot the r2
figure
bar(combination_results)
ylabel('R squared')
set(gca,'XTick',1:num_comb,'XTickLabels',...
    label_cell,'XTickLabelRotation',45,'TickLabelInterpreter','None')

% plot the rmse
figure
bar(combination_rmse)
ylabel('RMSE')
set(gca,'XTick',1:num_comb,'XTickLabels',...
    label_cell,'XTickLabelRotation',45,'TickLabelInterpreter','None')
%% Run the regression with shuffles

% get the selection vector (based on the OS cells)
selection_vector = vertcat(str.OSIpref)>0.25;
% define the target variables
response = 'ORIpref';

targets = {'ang_inL23_Vt','ang_inL4_Vt'};

% define the shuffle number and vector
shuffle_number = 100;
% shuffle_vector = [1 1];

% run the regression
[r2,coeff,shuffle_result,shuffle_prc]= GLM_fitting(str,selection_vector,response,...
    targets,shuffle_number);
%% Plot the shuffle results

close all

% r2 plot
bar(horzcat(r2,shuffle_result(1,:)))
hold on
errorbar(2:length(targets)+1,shuffle_result(1,:),shuffle_result(2,:),'ko')
ylabel('R squared')
set(gca,'XTick',1:length(targets)+1,'XTickLabels',...
    horzcat({'Original'},targets),'XTickLabelRotation',45,'TickLabelInterpreter','None')

% plot the coefficients
figure
bar(vertcat(coeff',shuffle_result(2:end,:)'))
ylabel('Coeff weight')
set(gca,'XTick',1:length(targets),'XTickLabels',...
    horzcat({'Original'},targets),'XTickLabelRotation',45,'TickLabelInterpreter','None')