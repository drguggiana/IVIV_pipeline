%% Predict gOSI based on input maps

clearvars
close all
matlabrc
Paths
%% Load the relevant files

% load the main structure
main_path = stage5_fullIVIV_path;
str = load(find_newer_file(main_path));
str = str.str;
%% Set a classifier to predict gOSI
close all
% get the number of cells
cell_num = length(str);
% load the gOSI
label_vector = [str.OSIpref];
% remove nans
nan_vector = isnan(label_vector);
label_vector = label_vector(~nan_vector);
% define the number of bins
n_bins = 4;
% bin
[~,~,label_vector] = histcounts(label_vector,n_bins,'BinMethod','fd');
% define the number of reps
rep_num = 10;
% allocate memory for the results
rep_cell = cell(rep_num,2);

% for the map types
for maptype = 1
    
    switch maptype
        case 1
            mtype = 'subpixel_excMap';
        case 2
            mtype = 'subpixel_inhMap';
    end
    
    % load the maps
    maps = cat(3,str.(mtype));
    % reshape
    maps = reshape(permute(maps,[3 1 2]),cell_num,[]);
    % remove nans
    maps = maps(~nan_vector,:);
    % normalize
    maps = normr_2(maps,1);
    % for all the reps
    for reps = 1:rep_num
        % set up the classifier
        %    model = fitrsvm(maps,label_vector,'Standardize',1,'KFold',5);
        model = fitcecoc(maps,label_vector,'KFold',5,'coding','onevsall','Prior','uniform');
        % predict
        y_fit = kfoldPredict(model);
        
        C = confusionmat(label_vector,y_fit);
        
        % store
        rep_cell{reps,maptype} = C;
        % plot the real and fit
        %         figure
        %         imagesc(C)
        %    scatter(label_vector,y_fit)
        %    hold on
        %    plot(y_fit)
    end
    % get the average confusion matrix
    average_C = mean(cat(3,rep_cell{:,maptype}),3);
    % plot the average confusion matrix
    figure
    imagesc(average_C)
    title(num2str(sum(diag(average_C))/sum(average_C(:))))
    axis square
%     set(gca
end

