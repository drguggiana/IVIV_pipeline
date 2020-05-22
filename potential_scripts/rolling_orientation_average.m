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
response = 'dir';

% define whether to plot the histogram
plot_histo = 1;

% define whether to run the shuffles
run_shuffle = 0;

% get the orientation and the parameter of interest
switch response
    case 'ori'
        orientation_vector = cat(1,str.ORIpref);
        selection_vector = cat(1,str.OSIpref)>0.25;
        axis_label = 'Orientation';
        plot_limit = 180;
        % define the window width (in degrees)
        window = 45;
    case 'dir'
        orientation_vector = cat(1,str.DIRpref);
        selection_vector = cat(1,str.DSIpref)>0.25;
        plot_limit = 360;
        axis_label = 'Direction';
        window = 65;
end
% get the selected traces
orientation_vector = orientation_vector(selection_vector);

% define the list of plots
plot_list = {'ang_exL23_Cx,ang_inL23_Cx','ang_exL23_Cy,ang_inL23_Cy',...
    'ang_exL23_Sx','ang_inL23_Sy','ang_exL23_Vt,ang_inL23_Vt','ang_exL23_Al,ang_inL23_Al',...
    'ang_exL23_Dx,ang_inL23_Dx','ang_exL23_Dy,ang_inL23_Dy'};
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
    % for maptype
    for mapnumber = 1:num_parameters
        % get the parameter for the rolling average
        plot_vector = parameter_vector(:,mapnumber);

        % run the function
        [rolling_orientation,rolling_ori_error,rolling_counts] = ...
            rolling_circular_average(orientation_vector,plot_vector,window,response);
        %% Generate a surrogate computation for a CI
        if run_shuffle == 1
            % define the number of shuffles
            shuffle_number = 500;
            % allocate memory to store the results
            shuffle_ori = zeros(shuffle_number,plot_limit);
            % take only the non-nan orientation and their values
            nonnan_ori = orientation_vector(~isnan(orientation_vector));
            nonnan_param = parameter_vector(~isnan(orientation_vector));

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

        %close all
        if plot_histo == 1
            subplot(2,1,1)
        end
        set(gcf, 'Position', [800, 200, 400, 300])
        set(gcf,'color','w');
        shadedErrorBar(1:plot_limit,rolling_orientation,rolling_ori_error,'transparent',1,...
            'lineprops',{'Color',cmap(mapnumber,:)})
        hold on
        %if the shuffle was activated
        if run_shuffle == 1
            shadedErrorBar(1:plot_limit,mean_shuffle,CI_shuffle,'transparent',1,'lineprops','k')
        end
        xlabel(strcat(axis_label,'(deg)'))
        ylabel('Parameter')
        title(strjoin({'Rolling',axis_label,'average,','window:',...
            num2str(window)},' '),'Interpreter','None')
        legend(strsplit(plot_list{plots},','),'Interpreter','None','Location','best')
        axis tight
    end
    if plot_histo == 1
        % also add a histogram
        subplot(2,1,2)
        bar(rolling_counts)
        set(gca,'XTick',20:20:plot_limit)
        xlabel('Orientation (deg)')
        ylabel('Nr points')
    end
end
autoArrangeFigures