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
%% Run the regression

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
set(gca,'XTick',1:length(targets)+1,'XTickLabels',...
    horzcat({'Original'},targets),'XTickLabelRotation',45,'TickLabelInterpreter','None')
%% OFF 
% % allocate memory for the actual data
% target_data = zeros(sum(selection_vector),length(targets));
% % for all the targets
% for target = 1:length(targets)
%     if contains(targets{target},'ang_')
%         
%         % parse the argument
%         field_name = targets{target}(1:end-3);
%         column_name = targets{target}(end-1:end);
%         % load the field
%         temp_var = vertcat(str(selection_vector).(field_name));
%         switch column_name
%             case 'Cx'
%                 target_data(:,target) = abs(temp_var(:,3)-temp_var(:,1));
%             case 'Cy'
%                 target_data(:,target) = abs(temp_var(:,4)-temp_var(:,2));
%             case 'Vt'
%                 target_data(:,target) = temp_var(:,8);
%             case 'Al'
%                 target_data(:,target) = temp_var(:,5);
%         end
%     elseif contains(targets{target},'frac_vert')
%         % parse the argument
%         field_name = targets{target}(1:9);
%         column_name = targets{target}(11:end);
%         
%         temp_var = vertcat(str(selection_vector).frac_vert);
%         switch column_name
%             case 'exL23'
%                 target_data(:,target) = sum(temp_var(:,3:5),2);
%             case 'exL4'
%                 target_data(:,target) = sum(temp_var(:,6:7),2);
%             case 'inL23'
%                 target_data(:,target) = sum(temp_var(:,19:21),2);
%             case 'inL4'
%                 target_data(:,target) = sum(temp_var(:,22:23),2);
%         end
%     elseif contains(targets{target},'custom')
%             temp_var = vertcat(str(selection_vector).ang_inL23);
%             x = abs(temp_var(:,3)-temp_var(:,1));
%             y = abs(temp_var(:,4)-temp_var(:,2));
%             target_data = sqrt(x.^2 + y.^2);
%     else
%             target_data(:,target) = vertcat(str(selection_vector).(targets{target}));
%             target_data(:,target) = target_data(randperm(sum(selection_vector)),target);
%     end
% end
% 
% % center the data
% target_data = (target_data-nanmean(target_data,1))./nanstd(target_data,0,1);
% 
% % load the predictor
% Y = vertcat(str(selection_vector).(dependent));
% Y = (Y-mean(Y))/std(Y);
% 
% % figure
% % 
% % corr(Y,target_data(:,1))
% % hold on
% % scatter(Y,target_data(:,2))
% % run the regression
% % mdl = fitrsvm(target_data,Y,'KernelFunction','gaussian');
% % mdl = fitglm(target_data,Y,[2 0 0;0 2 0],'Intercept',false); 
% mdl = fitglm(target_data,Y);
%% Plot the fit results
close all

figure
subplot(1,2,1)
scatter(Y,mdl.Fitted.Response)
title(mdl.Rsquared.Adjusted)
subplot(1,2,2)
bar(mdl.Coefficients.Estimate)
set(gca,'XTick',1:length(targets)+1,'XTickLabels',...
    horzcat({'Intercept'},targets),'XTickLabelRotation',45)
set(gca,'TickLabelInterpreter','none')

