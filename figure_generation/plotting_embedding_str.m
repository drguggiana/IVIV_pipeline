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
    switch plot_list{plots}
        case 'PCs'
            title_list = {'PC1_exc','PC2_exc','PC3_exc','PC1_inh','PC2_inh','PC3_inh',};
        case {'ang_exL23','ang_exL4','ang_inL23','ang_inL4'}
            data_in = cat(2,data_in(:,3)-data_in(:,1),data_in(:,5));
            title_list = {'centroidX_','alpha_'};
            title_list = cellfun(@strcat,title_list,...
                num2cell(repmat(plot_list(plots),1,size(title_list,2))),...
                'UniformOutput',false);
        case {'somaCenter','subpixel_soma'}
            title_list = {'soma x','soma y'};
        case 'cellID'
            title_list = {'setup'};
            data_in = cat(1,ones(47,1),zeros(100,1));
        case 'frac_vert'
            title_list = {'L4 exfrac','L4 infrac'};
            data_in = cat(1,str.frac_vert);
            data_in = [mean(data_in(:,6:7),2), mean(data_in(:,22:23),2)];
        case 'noise'
            data_in = log(abs(data_in));
            title_list = {'noise'};
        case 'pci'
            data_in = log(abs(data_in));
            title_list = plot_list(plots);
        otherwise
            title_list = plot_list(plots);
    end
    
    % define loop based on the dimensionality of data_in
    loop = size(data_in,2);
    
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