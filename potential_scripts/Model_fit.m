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
% reset the rng
rng(1)
%% Calculate the fit with different parameters

% get the selection vector (based on the OS cells)
selection_vector = vertcat(str.OSIpref)>0.25;
% define the target variables
response = 'ORIpref';

% targets = {'ang_inL23_Vt','ang_inL4_Cx','pialD','frac_vert_exL4','frac_vert_exL23',...
%     'ang_exL23_Vt','ang_inL4_Vt','frac_vert_inL4','frac_vert_inL23','nw_inL23_Vt'};
% targets = {'ang_inL23_Vt','max_ex_v_L23','max_ex_x_L23','max_ex_y_L23',...
%     'max_ex_v_L4','max_ex_x_L4','max_ex_y_L4','max_ex_v_L4','max_ex_x_L5','max_ex_y_L5'};
% targets = {'ang_inL23_Vt','max_in_v_L23','max_in_x_L23','max_in_y_L23',...
%     'max_in_v_L4','max_in_x_L4','max_in_y_L4','max_in_v_L4','max_in_x_L5','max_in_y_L5'};

targets = {'ang_inL23_Vt','ang_inL23_Cx','ang_inL23_Cy','ang_inL23_Sy','ang_inL23_Al','ang_inL23_Sx'};

% targets = {'ang_inL23_Vt','nw_inL23_Vt'};

% define the kernel function to use
kernel_function = 'rbf';

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
combination_loss = zeros(num_comb,1);
% allocate memory for the labels
label_cell = cell(num_comb,1);
% save the models
model_cell = cell(num_comb,2);
% for all the combinations
for combos = 1:num_comb
    % indicate the current combination
    fprintf(strjoin({'Current combination:',num2str(combos),...
        'out of',num2str(num_comb),'\r\n'},' '))
    % define the sub targets based on the combo
    % get the indexing vector
    idx_vector = comb(combos,:);
    idx_vector = idx_vector(idx_vector>0);
    sub_targets = targets(idx_vector);
    % run the regression
    [combination_results(combos),combination_coefficients{combos},...
        ~,~,combination_loss(combos),model_cell{combos,1},model_cell{combos,2}] =...
        SVR_fitting(str,selection_vector,response,...
        sub_targets,0,kernel_function);
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

% plot the loss
figure
bar(combination_loss)
ylabel('Loss')
set(gca,'XTick',1:num_comb,'XTickLabels',...
    label_cell,'XTickLabelRotation',45,'TickLabelInterpreter','None')

autoArrangeFigures
%% Plot a desired fit

close all
% define the target fot
% target_label = 'ang_inL23_Vt,frac_vert_inL4';
target_label = 'ang_inL23_Vt,ang_inL23_Sx';


% find the target model
target_model = model_cell{contains(label_cell,target_label),1};
target_nan_vector = model_cell{contains(label_cell,target_label),2};
target_coeff = combination_coefficients{contains(label_cell,target_label)};
% pull the data
real_data = vertcat(str(selection_vector).(response));
% filter the nandata with the vector
real_data = real_data(target_nan_vector);
% real_data = (real_data-mean(real_data))./std(real_data);

% get the target model data
model_data = kfoldPredict(target_model);

% add the mean and std info from the real data
% model_data = model_data.*std(real_data) + mean(real_data);

% calculate a simple linear fit
lin_mdl = fitlm(model_data,real_data);

figure
scatter(real_data,model_data,'k','filled')
hold on
plot(real_data,predict(lin_mdl,real_data))
xlabel('Real data')
ylabel('Predicted data')
title(strjoin({'Model: ',target_label,'r2:',num2str(lin_mdl.Rsquared.Adjusted)},' '),'Interpreter','None')
axis equal

if ~isempty(target_coeff)
    figure
    bar(target_coeff)
    targets = strsplit(target_label,',');
    ylabel('Coeff weight')
    set(gca,'XTick',1:length(targets),'XTickLabels',...
        targets,'XTickLabelRotation',45,'TickLabelInterpreter','None')
end
autoArrangeFigures
%% Run the regression with shuffles

% get the selection vector (based on the OS cells)
selection_vector = vertcat(str.OSIpref)>0.25;
% define the target variables
response = 'ORIpref';

targets = {'ang_inL23_Vt','ang_inL23_Sx'};

% define the shuffle number and vector
shuffle_number = 100;
% shuffle_vector = [1 1];

% run the regression
[r2,coeff,shuffle_result,shuffle_prc]= SVR_fitting(str,selection_vector,response,...
    targets,shuffle_number,kernel_function);
%% Plot the shuffle results

close all

% r2 plot
bar(horzcat(r2,shuffle_result(1,:)))
hold on
errorbar(2:length(targets)+1,shuffle_result(1,:),shuffle_prc(1,:),'ko')
ylabel('R squared')
set(gca,'XTick',1:length(targets)+1,'XTickLabels',...
    horzcat({'Original'},targets),'XTickLabelRotation',45,'TickLabelInterpreter','None')

if ~isempty(target_coeff)
    % plot the coefficients
    figure
    bar(vertcat(coeff',shuffle_result(2:end,:)'))
    ylabel('Coeff weight')
    set(gca,'XTick',1:length(targets)+1,'XTickLabels',...
        horzcat({'Original'},targets),'XTickLabelRotation',45,'TickLabelInterpreter','None')
end
autoArrangeFigures