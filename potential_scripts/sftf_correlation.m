%% Correlation schemes
%% load the paths and clean up
clearvars
close all

Paths

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Check different correlation functions

sf = vertcat(str.SF);
tf = vertcat(str.TF);

% define the target parameter to predict
target_parameter = tf;
% get a vector with the non-nan data
not_nan = ~isnan(target_parameter);

% filter the data
target_parameter = target_parameter(not_nan);

% try a classifier
[~,~,label_vector] = unique(target_parameter);
% label_vector(label_vector==3) = nan;

% select a parameter
% plot_list = {'ang_exL23_Cx','ang_inL23_Cx','ang_exL23_Cy','ang_inL23_Cy',...
%     'ang_exL23_Sx','ang_inL23_Sy','ang_exL23_Vt','ang_inL23_Vt','ang_exL23_Al','ang_inL23_Al',...
%     'ang_exL23_Dx','ang_inL23_Dx','ang_exL23_Dy','ang_inL23_Dy'};
% plot_list = {'ang_exL4_Cx','ang_inL4_Cx','ang_exL4_Cy','ang_inL4_Cy',...
%     'ang_exL4_Vt','ang_inL4_Vt','ang_exL4_Al','ang_inL4_Al',...
%     'ang_exL4_Rx','ang_inL4_Rx','ang_exL4_Ry','ang_inL4_Ry'};
% plot_list = {'ang_exL5_Cx','ang_inL5_Cx','ang_exL5_Cy','ang_inL5_Cy',...
%     'ang_exL5_Vt','ang_inL5_Vt','ang_exL5_Al','ang_inL5_Al',...
%     'ang_exL5_Rx','ang_inL5_Rx','ang_exL5_Ry','ang_inL5_Ry'};
% plot_list = {'span_L23_ex','span_L23_in','span_L4_ex','span_L4_in','span_L5_ex','span_L5_in'};
% plot_list = {'frac_vert_exL23','frac_vert_inL23','frac_vert_exL4','frac_vert_inL4',...
%     'frac_vert_exL5','frac_vert_inL5'};
% plot_list = {'max_ex_v_L23,max_in_v_L23','max_ex_x_L23,max_in_x_L23','max_ex_y_L23,max_in_y_L23',...
%     'max_ex_v_L4,max_in_v_L4','max_ex_x_L4,max_in_x_L4','max_ex_y_L4,max_in_y_L4',...
%     'max_ex_v_L5,max_in_v_L5','max_ex_x_L5,max_in_x_L5','max_ex_y_L5,max_in_y_L5'};

plot_list = {'frac_vert_exL23','frac_vert_inL23','frac_vert_exL4','frac_vert_inL4',...
    'frac_vert_exL5','frac_vert_inL5','ang_exL23_Cx','ang_inL23_Cx','ang_exL23_Cy','ang_inL23_Cy',...
    'span_L23_ex','span_L23_in','span_L4_ex','span_L4_in','span_L5_ex','span_L5_in',...
    'pialD'};



% just select all in vivo
selection_vector = not_nan;


% define the number of repeats
repeat_number = 100;

% allocate memory for the performances
performances = zeros(length(plot_list),repeat_number);

% allocate memory for the shuffle results
shuffles = zeros(length(plot_list),repeat_number);


% for all the items in plot_list
for plots = 1:length(plot_list)
    
    fprintf(strjoin({'Current plot:',num2str(plots),...
        'out of',num2str(length(plot_list)),'\r\n'},' '))
    parameter_vector = variable_selector(plot_list{plots},str,selection_vector);
    
    % normalize the parameter
    parameter_vector = normr_2(parameter_vector,0);
    
    [parameter_vector,y] = nan_remover(parameter_vector,label_vector);
    
    % for all the repeats
    for repeat = 1:repeat_number
        Mdl = fitcecoc(parameter_vector,y,'KFold',5,'Learners','naivebayes');
        
        y_pred = kfoldPredict(Mdl);
        
        performances(plots,repeat) = sum(y_pred==y)./length(y);
        
        % shuffle and calculate performance
        shuffle_parameter = parameter_vector(randperm(length(parameter_vector)));
        Mdl_shuffle = fitcecoc(shuffle_parameter,y,'KFold',5,'Learners','naivebayes');
        y_pred_shuffle = kfoldPredict(Mdl_shuffle);
        shuffles(plots,repeat) = sum(y_pred_shuffle==y)./length(y);
    end
end
%% Plot the results

close all

% calculate the average and std shuffle and performance
mean_performance = mean(performances,2);
sem_performance = std(performances,0,2);

mean_shuffle = mean(shuffles,2);
sem_shuffle = std(shuffles,0,2);

figure

% get the x range
x_range = 1:length(mean_performance);

bar(mean_performance)
hold on
errorbar(x_range,mean_performance,sem_performance,'.b')

errorbar(x_range,mean_shuffle,sem_shuffle,'k.')

set(gca,'XTick',x_range,'XTickLabels',plot_list,...
    'XTickLabelRotation',45,'TickLabelInterpreter','None')