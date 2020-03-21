function plotting_embedding2(reduced_data, str, plot_selector, non_nan_cells,  correlation_values,dim, varargin)

% for all the elements in plot selector
for el = 1:length(plot_selector)
    % get the number of cells (in the reduced data)
    cell_num = size(reduced_data,1);
    
    switch plot_selector(el)
        case 1
            % pial depth
            parameter = round([str(non_nan_cells).pialD]);
            data_title = 'Pial Depth';
        case 2
            % morphological parameters
            iv_param = {str(non_nan_cells).morph};
            data_title = 'Morpho parameter';
        case 3
            % orientation preference
            iv_param = {str(non_nan_cells).Oripref}';
            data_title = 'Orientation pref';
        case 4
            % direction preference
            iv_param = {str(non_nan_cells).Dirpref}';
            data_title = 'Direction pref';
        case 5
            % spontaneous activity
            iv_param = {str(non_nan_cells).iv_spon}';
            data_title = 'Spont act';
        case 6
            % population coupling
            iv_param = {str(non_nan_cells).iv_popcop}';
            data_title = 'Pop coupling';
        case 7
            % morpho basal | ex_input correlation
            iv_param = correlation_values(:,1);
            data_title = 'Ex_input to apical corr';
        case 8
            % morpho basal | in_input correlation
            iv_param = correlation_values(:,2);
            data_title = 'Inh_input to apical corr';
        case 9
            % setup COMPLETELY ARBITRARY TO SIMON'S DATA SO WATCH OUT
            parameter = [zeros(1,45),ones(1,112)];
            data_title = 'Setup';
        case 10
            % time
%             parameter = 1:157;
            cell_date = cellfun(@(x) str2double(x(3:end)),{str.cellName});
            [~, parameter] = sort(cell_date);
            data_title = 'Time';
        case 11
            % hemisphere
            parameter = [str(:).hemisphere];
            data_title = 'Hemisphere';
        case 12
            % slice orientation
            parameter = [str(:).sliceOri];
            data_title = 'Slice orientation';
         case 13
             %ODI
            iv_param = {str(non_nan_cells).ODI}';
             data_title = 'ODI';
         case 14
             %OSI
            iv_param = {str(non_nan_cells).OSI}';
             data_title = 'OSI';
         case 15
             %DSI
            iv_param = {str(non_nan_cells).DSI}';
             data_title = 'DSI';
        case 16
            %Peak Ca
            iv_param = cellfun(@log,{str(non_nan_cells).Peak_Ca}', 'UniformOutput',0);
            data_title = 'Peak Calcium';
        case 17
            %SF
            iv_param = {str(non_nan_cells).SF}';
            data_title = 'SF';
        case 18
            %TF
            iv_param = {str(non_nan_cells).TF}';
            data_title = 'TF';
        case 28
            %opp Resp pref
            iv_param = {str(non_nan_cells).oppResp}';
            data_title = 'Opp Resp Ori';
        case 29
            %opp Resp pref
            iv_param = {str(non_nan_cells).sigma_tuning}';
            data_title = 'Sigma tuning curve';
            
        otherwise
            % arbitrary parameter
            parameter = varargin{1}(:,plot_selector(el)-29);
            data_title = varargin{2}{plot_selector(el)-29};
            % % input map clusters
            % parameter = idx_input;
            
            % % TMD clusters
            % parameter = clusters;
             
    end
    % if the size of the parameter vector is higher than the reduced data
    if size(reduced_data,1)<size(parameter,1)
        % remove the non nan cells
        parameter = parameter(non_nan_cells);
    end
    
    % process the invivo or morpho parameters
    if any(plot_selector(el)==[2:8 13:29])
        parameter = zeros(cell_num,1);
        for cells = 1:cell_num
            if ~isempty(iv_param{cells}) == 1
                parameter(cells) = round(iv_param{cells}(1).*100);
               %CHANGED SW 11/08/2019
               % parameter(cells) = round(nanmean(iv_param{cells}).*100)
               
               %parameter(cells) = round(nanmean(iv_param{cells}))
            else
                parameter(cells) = NaN;
            end
        end
    end
    % define the color scale
   
    color_number = max(parameter)-min(parameter)+1;
    color_indexer = parameter-min(parameter)+1;
    cmap = parula(color_number);

    figure
    set(gcf,'color','w');
    % for all the cells
    for cells = 1:cell_num
        if dim == 2
            if ~isnan(parameter(cells))
                plot(reduced_data(cells,1),reduced_data(cells,2),'o','MarkerFaceColor',cmap(color_indexer(cells),:),'MarkerEdgeColor',cmap(color_indexer(cells),:))
                hold on
            else
                %plot(reduced_data(cells,1),reduced_data(cells,2),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
               plot(reduced_data(cells,1),reduced_data(cells,2),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.8 0.8 0.8]);
                hold on
            end
        else
            if ~isnan(parameter(cells))
                scatter3(reduced_data(cells,1),reduced_data(cells,2),reduced_data(cells,3),'o','MarkerFaceColor',cmap(color_indexer(cells),:),'MarkerEdgeColor',cmap(color_indexer(cells),:))
                hold on
            else
                %plot(reduced_data(cells,1),reduced_data(cells,2),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
               scatter3(reduced_data(cells,1),reduced_data(cells,2),reduced_data(cells,3),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.8 0.8 0.8]);
                hold on
            end
        end
    end
    title(data_title, 'interpreter', 'none')
    axis square
end

autoArrangeFigures