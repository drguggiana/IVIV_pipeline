function [r2,coeff,varargout]= GLM_fitting(str,selection_vector,dependent,targets,varargin)

% if there is a shuffle vector
if length(varargin) >= 1
    shuffle_number = varargin{1};
else
    shuffle_number = 0;
end

if length(varargin) >= 2
    cross_validate = varargin{2};
else
    cross_validate = false;
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
    else
            target_data(:,target) = vertcat(str(selection_vector).(targets{target}));
            target_data(:,target) = target_data(randperm(sum(selection_vector)),target);
    end
    
end

% load the dependent
Y = vertcat(str(selection_vector).(dependent));
Y = (Y-mean(Y))/std(Y);

% remove missing values
nan_vector = sum(isnan(target_data),2)==0;
target_data = target_data(nan_vector,:);
Y = Y(nan_vector);

% center the data
target_data = (target_data-nanmean(target_data,1))./nanstd(target_data,0,1);

% if shuffles are desired , shuffle according to the shuffling vector
if shuffle_number > 0
    % allocate memory to store the shuffles
    shuffle_out = zeros(shuffle_number,size(target_data,2)+1,size(target_data,2));
    % for all the shuffles
    for shuff = 1:shuffle_number
%         % get the indexes for the columns to be shuffled
%         shuff_cols_number = sum(shuffle_vector);
        % allocate the permutation matrix
        temp_data = target_data;
        % for all the shuffle columns
        for col = 1:size(target_data,2)
            % generate the permutations
            temp_data(:,col) = temp_data(randperm(size(temp_data,1)),col);
            % perform the fit
            mdl = fitglm(temp_data,Y,'linear','Distribution','normal','intercept',false);
            % save the coefficients and fit
            shuffle_out(shuff,1,col) = mdl.Rsquared.Adjusted;
            shuffle_out(shuff,2:end,col) = mdl.Coefficients.Estimate;
        end

    end
    % output the results
    shuffle_result = mean(shuffle_out,1);
    shuffle_prc = prctile(shuffle_out,95,1)-shuffle_result;
    varargout{1} = squeeze(shuffle_result);
    varargout{2} = squeeze(shuffle_prc);
end

% if cross val is on, cross validate
if cross_validate
    fcn = @(Xtr, Ytr, Xte) predict(...
    GeneralizedLinearModel.fit(Xtr,Ytr,'linear','distr','normal'), ...
    Xte);

    % perform cross-validation, and return average MSE across folds
    varargout{3} = sqrt(crossval('mse', target_data, Y, 'Predfun',fcn, 'kfold',10));
    
end
% calculate the model and output
mdl = fitglm(target_data,Y,'linear','Distribution','normal','intercept',false);
r2 = mdl.Rsquared.Adjusted;
coeff = mdl.Coefficients.Estimate;

varargout{4} = mdl;
<<<<<<< HEAD
varargout{5} = nan_vector;
=======
varargout{5} = nan_vector;

>>>>>>> d978a27a5ad3f9256da7db2ff072b009dafaa506
