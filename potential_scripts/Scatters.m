%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Calculate the rolling average for orientation
close all

% define whether to do orientation or direction
response = 'ori';

% define whether to run the shuffles
run_shuffle = 0;

% get the orientation and the parameter of interest
switch response
    case 'ori'
        angle_vector = cat(1,str.ORIpref);
%         selection_vector = cat(1,str.OSIpref)>0.25;
        selection_vector = ones(147,1)==1;
%         selection_vector = cat(1,str.DSIpref)>0.25;

        axis_label = 'Orientation';
        plot_limit = 180;
        % define the window width (in degrees)
        window = 45;
    case 'dir'
        angle_vector = cat(1,str.DIRpref);
        selection_vector = cat(1,str.DSIpref)>0.25;
        plot_limit = 360;
        axis_label = 'Direction';
        window = 65;
end
% get the selected traces
% angle_vector = angle_vector(selection_vector);
angle_vector = selection_vector;

% define the list of plots
% plot_list = {'max_ex_x_L23,ang_exL23_Cx','max_in_x_L23,ang_inL23_Cx','ang_exL23_Cy,ang_inL23_Cy',...
%     'ang_exL23_Vt,ang_inL23_Vt','ang_exL23_Al,ang_inL23_Al',...
%     'ang_exL23_Dx,ang_inL23_Dx','ang_exL23_Dy,ang_inL23_Dy'};

plot_list = {'OSIpref,DSIpref'};

% plot_list = {'max_in_x_L23,frac_vert_inL4','max_ex_x_L23,frac_vert_exL4',...
%     'max_ex_y_L23,frac_vert_exL5','max_in_y_L23,frac_vert_inL5',...
%     'max_ex_v_L4,ang_exL23_Dx','max_in_v_L4,ang_inL23_Dx'};
% plot_list = {'ang_exL4_Cx,ang_inL4_Cx','ang_exL4_Cy,ang_inL4_Cy',...
%     'ang_exL4_Vt,ang_inL4_Vt','ang_exL4_Al,ang_inL4_Al',...
%     'ang_exL4_Rx,ang_inL4_Rx','ang_exL4_Ry,ang_inL4_Ry'};
% plot_list = {'ang_exL5_Cx,ang_inL5_Cx','ang_exL5_Cy,ang_inL5_Cy',...
%     'ang_exL5_Vt,ang_inL5_Vt','ang_exL5_Al,ang_inL5_Al',...
%     'ang_exL5_Rx,ang_inL5_Rx','ang_exL5_Ry,ang_inL5_Ry'};
% plot_list = {'span_L23_ex,span_L23_in','span_L4_ex,span_L4_in','span_L5_ex,span_L5_in'};
% plot_list = {'frac_vert_exL23,frac_vert_inL23','frac_vert_exL4,frac_vert_inL4',...
%     'frac_vert_exL5,frac_vert_inL5'};
% plot_list = {'max_ex_v_L23,max_in_v_L23','max_ex_x_L23,max_in_x_L23','max_ex_y_L23,max_in_y_L23',...
%     'max_ex_v_L4,max_in_v_L4','max_ex_x_L4,max_in_x_L4','max_ex_y_L4,max_in_y_L4',...
%     'max_ex_v_L5,max_in_v_L5','max_ex_x_L5,max_in_x_L5','max_ex_y_L5,max_in_y_L5'};

% get the number of plots
plot_number = numel(plot_list);

% for all the plots
for plots = 1:plot_number
    % get the parameter vector
    parameter_vector = variable_selector(plot_list{plots},str,selection_vector);
    % get the number of parameters
    num_parameters = size(parameter_vector,2);
    % get a colormap
    cmap = lines(num_parameters);


    figure
%     % for maptype
%     for mapnumber = 1:num_parameters
        % get the parameter for the rolling average
        plot_x = parameter_vector(:,1);
        plot_y = parameter_vector(:,2);

        
        %% Generate a surrogate computation for a CI
        if run_shuffle == 1
            % define the number of shuffles
            shuffle_number = 500;
            % allocate memory to store the results
            shuffle_ori = zeros(shuffle_number,plot_limit);
            % take only the non-nan orientation and their values
            nonnan_ori = angle_vector(~isnan(angle_vector));
            nonnan_param = parameter_vector(~isnan(angle_vector));

            % for all the shuffles
            for shuffles = 1:shuffle_number
                % randomize the parameter_vector
                rand_param = nonnan_param(randperm(length(nonnan_param)));
                % run the function
                [shuffle_ori(shuffles,:),~] = rolling_circular_average(nonnan_ori,rand_param,window,response);
            end
            % get the mean and bounds
            mean_shuffle = nanmean(shuffle_ori,1);
            CI_shuffle = cat(1,abs(prctile(shuffle_ori,5,1)-mean_shuffle),prctile(shuffle_ori,95,1)-mean_shuffle);
        end
        %% Plot the results


        set(gcf, 'Position', [800, 200, 400, 300])
        set(gcf,'color','w');
        scatter(plot_x,plot_y,30,angle_vector)
        hold on
        %if the shuffle was activated
        if run_shuffle == 1
            shadedErrorBar(1:plot_limit,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
        end
        labels = strsplit(plot_list{plots},',');
        xlabel(labels{1},'Interpreter','None')
        ylabel(labels{2},'Interpreter','None')
        [plot_x,plot_y] = nan_remover(plot_x,plot_y);
        [rho,pval] = corr(plot_x,plot_y);
        title(strjoin({num2str(rho),num2str(pval)},' '),'Interpreter','None')
        axis equal
        axis tight
%     end

end
autoArrangeFigures