function [r2,coeff,varargout]= SVR_fitting(str,selection_vector,dependent,targets,varargin)

% if there is a shuffle vector
if length(varargin) >= 1
    shuffle_number = varargin{1};
else
    shuffle_number = 0;
end

if length(varargin) >= 2
    kernel_function = varargin{2};
else
    kernel_function = 'linear';
end


% allocate memory for the actual data
target_data = zeros(sum(selection_vector),length(targets));
% for all the targets
for target = 1:length(targets)
    if contains(targets{target},'ang_')
        
        % parse the argument
        field_name = targets{target}(1:end-3);
        column_name = targets{target}(end-1:end);
        % load the field
        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'Cx'
                target_data(:,target) = abs(temp_var(:,3)-temp_var(:,1));
            case 'Cy'
                target_data(:,target) = abs(temp_var(:,4)-temp_var(:,2));
            case 'Vt'
                target_data(:,target) = temp_var(:,8);
            case 'Al'
                target_data(:,target) = temp_var(:,5);
        end
    elseif contains(targets{target},'nw_')
        % parse the argument
        field_name = strcat('ang',targets{target}(3:end-3),'_nonweighted');
        column_name = targets{target}(end-1:end);
        % load the field
        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'Cx'
                target_data(:,target) = abs(temp_var(:,3)-temp_var(:,1));
            case 'Cy'
                target_data(:,target) = abs(temp_var(:,4)-temp_var(:,2));
            case 'Vt'
                target_data(:,target) = temp_var(:,8);
            case 'Al'
                target_data(:,target) = temp_var(:,5);
        end
    elseif contains(targets{target},{'frac_vert','frac_horz'})
        % parse the argument
        field_name = targets{target}(1:9);
        column_name = targets{target}(11:end);
        
        temp_var = vertcat(str(selection_vector).(field_name));
        switch column_name
            case 'exL23'
                target_data(:,target) = sum(temp_var(:,3:5),2);
            case 'exL4'
                target_data(:,target) = sum(temp_var(:,6:7),2);
            case 'inL23'
                target_data(:,target) = sum(temp_var(:,19:21),2);
            case 'inL4'
                target_data(:,target) = sum(temp_var(:,22:23),2);
        end
    elseif contains(targets{target},'custom')
        temp_var = vertcat(str(selection_vector).ang_inL23);
        x = abs(temp_var(:,3)-temp_var(:,1));
        y = abs(temp_var(:,4)-temp_var(:,2));
        target_data = sqrt(x.^2 + y.^2);
    elseif contains(targets{target},'max_')
        % parse the argument
        field_name = targets{target}(1:6);
        coord_name = targets{target}(8);
        layer_name = targets{target}(10:end);
        
        % get the data
        temp_var = reshape(vertcat(str(selection_vector).(field_name)),[],3,3);
        switch coord_name
            case 'v'
                coord = 1;
            case 'x'
                coord = 2;
            case 'y'
                coord = 3;
        end
        switch layer_name
            case 'L23'
                layer = 1;
            case 'L4'
                layer = 2;
            case 'L5'
                layer = 3;
        end
        
        % store the data
        target_data(:,target) = squeeze(temp_var(:,coord,layer));
        
    else
        target_data(:,target) = vertcat(str(selection_vector).(targets{target}));
        target_data(:,target) = target_data(randperm(sum(selection_vector)),target);
    end
    
end

% load the dependent
Y = vertcat(str(selection_vector).(dependent));
% Y = (Y-mean(Y))/std(Y);

% % remove missing values
nan_vector = sum(isnan(target_data),2)==0;
target_data = target_data(nan_vector,:);
Y = Y(nan_vector);

% % center the data
% target_data = (target_data-nanmean(target_data,1))./nanstd(target_data,0,1);

% if shuffles are desired , shuffle according to the shuffling vector
if shuffle_number > 0
    % allocate memory to store the shuffles
    shuffle_out = zeros(shuffle_number,size(target_data,2)+1,size(target_data,2));
    % for all the shuffles
    for shuff = 1:shuffle_number
        
        % display the current shuffle
        fprintf(strjoin({'Current shuffle:',num2str(shuff),'out of',...
            num2str(shuffle_number),'\r\n'},' '))
%         % get the indexes for the columns to be shuffled
%         shuff_cols_number = sum(shuffle_vector);
        % allocate the permutation matrix
        temp_data = target_data;
        % for all the shuffle columns
        for col = 1:size(target_data,2)
            % generate the permutations
            temp_data(:,col) = temp_data(randperm(size(temp_data,1)),col);
            % perform the fit
%             mdl = fitglm(temp_data,Y,'linear','Distribution','normal','intercept',false);
            mdl = fitrsvm(temp_data,Y,'KFold',40,'KernelFunction',kernel_function,...
                'Standardize',true,'KernelScale',1);
%             mdl = fitrgp(temp_data,Y,'Standardize',true,'Leaveout','on');

            % save the coefficients and fit
            lm = fitlm(Y,kfoldPredict(mdl));
            temp_r2 = lm.Rsquared.Adjusted;
            shuffle_out(shuff,1,col) = temp_r2;
            
            % if the kernel function is not linear, skip
            if ~strcmp(kernel_function,'linear')
                continue
            end
            temp_coeff = cellfun(@(x) x.Beta,mdl.Trained,'UniformOutput',false);
            temp_coeff = mean(horzcat(temp_coeff{:}),2);
            shuffle_out(shuff,2:end,col) = temp_coeff;
        end

    end
    % output the results
    shuffle_result = mean(shuffle_out,1);
    shuffle_prc = prctile(shuffle_out,95,1)-shuffle_result;
    varargout{1} = squeeze(shuffle_result);
    varargout{2} = squeeze(shuffle_prc);
end

% calculate the actual fit
% mdl = fitrgp(target_data,Y,'Standardize',true,'Leaveout','on');
mdl = fitrsvm(target_data,Y,'KFold',40,'KernelFunction',kernel_function,...
    'Standardize',true,'KernelScale',1);
varargout{3} = kfoldLoss(mdl);

% if the kernel function is not linear, skip
if strcmp(kernel_function,'linear')
    coeff = cellfun(@(x) x.Beta,mdl.Trained,'UniformOutput',false);
    coeff = mean(horzcat(coeff{:}),2);
else
    coeff = [];
end
% calculate the rsquared
lm = fitlm(Y,kfoldPredict(mdl));
r2 = lm.Rsquared.Adjusted;

varargout{4} = mdl;
varargout{5} = nan_vector;
