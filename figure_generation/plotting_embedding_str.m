function plotting_embedding_str(reduced_data,str,plot_list,varargin)


% get the number of plots
plot_number = length(plot_list);
% get the log scale vector if available
if length(varargin) >= 1
    if plot_number > length(varargin{1})
        log_vector = ones(plot_number,1).*varargin{1};
    else
        log_vector = varargin{1};
    end
else
    log_vector = zeros(plot_number,1);
end
% get the normalization vector if available
if length(varargin) >= 2
    if plot_number > length(varargin{2})
        norm_vector = ones(plot_number,1).*varargin{2};
    else
        norm_vector = varargin{2};
    end
else
    norm_vector = zeros(plot_number,1);
end
% get the colormap if available
if length(varargin) >= 3
    cmap = varargin{3};
else
    cmap = 'parula';
end

% for all the elements in the plot_list
for plots = 1:plot_number
    
    % get the data from the selected field
    data_in = cat(1,str.(plot_list{plots}));
    
    % if it's the PCs, split and plot all
    if strcmp(plot_list{plots},'PCs') == 1
        % set the loop counter
        loop = size(data_in,2);
        title_list = {'PC1_exc','PC2_exc','PC3_exc','PC1_inh','PC2_inh','PC3_inh',};
    elseif contains(plot_list{plots},'ang_')
        loop = 1;
        data_in = data_in(:,5);
        title_list = plot_list(plots);
    else
        loop = 1;
        title_list = plot_list(plots);
    end
    
    % for all the loop values
    for count = 1:loop
        % if there's only one iteration, just concatenate
        if loop > 1
            data = data_in(:,count);
        else
            data = data_in;
        end

        % normalize if the flag indicates so
        if norm_vector(plots) == 1
            data = 1 + normr_2(data);
        end
        % if log scale, turn to log
        if log_vector(plots) == 1
            data = log(data);
        end

        % split between nan and not nan
        nan_data = reduced_data(isnan(data),:);
        value_data = reduced_data(~isnan(data),:);

        figure
        set(gcf,'color','w');
        scatter(nan_data(:,1),nan_data(:,2),30,[0.8, 0.8 0.8])
        hold on
        scatter(value_data(:,1),value_data(:,2),30,data(~isnan(data)),'filled')
        colormap(cmap)

        title(title_list{count}, 'interpreter', 'none')
        axis square
    end
end