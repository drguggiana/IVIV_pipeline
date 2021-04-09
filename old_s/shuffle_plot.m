function [mean_shuffle CI_shuffle ] = shuffle_plot(angle_vector,plot_vector)

shuffle_number = 500;
            % allocate memory to store the results
%             shuffle_ori = zeros(shuffle_number,plot_limit);
            shuffle_ori = zeros(shuffle_number,1);

            % take only the non-nan orientation and their values
%             nonnan_ori = angle_vector(~isnan(angle_vector));
%             nonnan_param = parameter_vector(~isnan(angle_vector));
            [nonnan_ori,nonnan_param] = nan_remover(angle_vector,plot_vector);

            % for all the shuffles
            for shuffles = 1:shuffle_number
                % randomize the parameter_vector
                rand_param = nonnan_param(randperm(length(nonnan_param)));
                % run the function
%                 [shuffle_ori(shuffles,:),~] = rolling_circular_average(nonnan_ori,rand_param,window,response);

                [sorted_angle,sort_idx] = sort(nonnan_ori);
                sorted_points = rand_param(sort_idx);
                shuffle_ori(shuffles) = fit_sine(sorted_angle',sorted_points',0);
            end
            % get the mean and bounds
            mean_shuffle = nanmean(shuffle_ori,1);
            CI_shuffle = cat(1,abs(prctile(shuffle_ori,5,1)-mean_shuffle),...
                prctile(shuffle_ori,95,1)-mean_shuffle);
            % store them for later use
%             shuffles_temp(mapnumber,1) = mean_shuffle;
%             shuffles_temp(mapnumber,2:3) = CI_shuffle;
end