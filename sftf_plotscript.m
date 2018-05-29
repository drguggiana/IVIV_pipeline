%% Load SF and TF data and plot as arrow plots

%% Clean up
clearvars
close all
%% Load the excel file

%define the path
xls_path = 'R:\Share\Simon\Drago_Volker_Simon\Invivo_all_animals_expID.xlsx';

%load the file
[~,~,xls_raw] = xlsread(xls_path);
%% Load the actual data

%define the target path
main_path = 'I:\Simon Weiler\AnalyzedData\';

%get a list of folders within it
folder_list = dir(main_path);
%get rid of the dot folders
folder_list = folder_list(3:end);
%also of the non folder items
folder_list = folder_list(cat(1,folder_list(:).isdir));

%modify the folder list based on the excel spreadsheet
[~,ia,ic] = intersect(xls_raw(2:end,1),{folder_list(:).name});
%get the ordered binary vector for including files
include_vec = cat(1,xls_raw{2:end,11});
include_vec = include_vec(ia)==1;
folder_list = folder_list(include_vec);

%get the number of folders
folder_num = length(folder_list);
%go folder by folder loading the data
%allocate memory for the actual data
sftf_cell = cell(folder_num,1);
%allocate memory to store the experiment name
name_cell = cell(folder_num,1);

%for all the folders
for folders = 1:folder_num
    %define the target path
    tar_path = strcat(folder_list(folders).folder,'\',folder_list(folders).name,...
        '\*\Analysis\sftf_analysis');
    %get the file name
    tar_name = dir(strcat(tar_path,'\*.mat'));
    %load the file
    sftf_cell{folders} = load(fullfile(tar_name(1).folder,tar_name(1).name));
    %save the experiment name
    name_cell{folders} = str2double(tar_name.name(10:14));
    
    
end
%get the actual sf, tf and dirs
sf_vec = sftf_cell{1}.ids.stim_SFs;
tf_vec = sftf_cell{1}.ids.stim_TFs;
dir_vec = sftf_cell{1}.ids.stim_dirs;
%% Format the data

%allocate memory for the sf tf dir matrices for each experiment
sftf_mat = cell(folder_num,2);
%elements are in dir, sf, tf order, assume correct order, otherwise will
%have to load sequence
%also store the PSTHs in a separate cell
psth_cell = cell(folder_num,1);

%define the threshold p value
p_val = 0.05;

%for all the folders
for folders = 1:folder_num
    %get the experiment name and save it
    sftf_mat{folders,1} = name_cell{folders};
    
    %get the corresponding matrix, normalize within the experiment to 8
    %bits
    sftf_mat{folders,2} = normr_2(permute(sftf_cell{folders}.peaks.average_trace_delta_integ,[1 3 4 2])).*255;
    %filter the matrix with the ANOVA information
    %get the ANOVA info for each SF and TF
%     sftf_anova = sftf_cell{folders}.peaks.ANOVA_p_singleSF_TF_only;    
%     %redimension the anova matrix to match the data matrix (add dir dim)
%     sftf_anova = permute(repmat(sftf_anova,1,1,1,length(dir_vec)),[1 4 2 3]);
    sftf_anova = sftf_cell{folders}.peaks.ANOVA_p_all_conditions;

    %filter the data matrix
%     sftf_mat{folders,2}(sftf_anova>p_val) = NaN;
    sftf_mat{folders,2}(sftf_anova>p_val,:,:,:) = NaN;

    
    %save the psth separately
    psth_cell{folders} = sftf_cell{folders}.peaks.mean_full_PSTH_data;
         
end

%resort the mat based on exp name
[~,idx] = sort(cat(1,sftf_mat{:,1}));
sftf_mat = sftf_mat(idx,:);
%concatenate all the cells together
sftf_all = cat(1,sftf_mat{:,2});
psth_all = cat(1,psth_cell{:});
%get the number of cells
cell_num = size(sftf_all,1);

%convert the cell to a structure
sftf_mat = cell2struct(sftf_mat,{'expName','data'},2);
%% Save the cell as a structure including the experiment ID

%define the save path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\sftf_out\';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_sftfData.mat');

%save the file
save(fullfile(save_path,save_name),'sftf_mat');
%% Run the plotting function for a subset of cells

close all
%define the number of cells to plot
tar_cells = 10;
%get the vector of maps to plot
map_vec = randperm(cell_num,tar_cells);
%for all the desired cells
for cells = map_vec
    %create a figure to plot things into
    h = figure;
    
    %run the plotting function
    sftf_plot(squeeze(sftf_all(cells,:,:,:)),sf_vec,tf_vec,dir_vec,h,squeeze(psth_all(cells,:,:,:)))
%     sftf_plot(squeeze(sftf_all(cells,:,:,:)),sf_vec,tf_vec,dir_vec,h)
end
%% Plot the iviv cells only

% close all
iviv={1:3;1;1;1:2;1:3;1;1:2;1:3;1:3;1:2;1:3;1:3;1:4;1:3;1:3;1:2;1:3;1:3;1:2;1:3;1:3;1:2;1:2;1:3;1:3;1;1:3;1:3;1;1:3;1:2;1:3};
% 
% %initialize subplot counter
% sub_count = 1;
% %create a figure to plot things into
% h = figure;
% %for all the folders
% for folders = 1:folder_num
%     %get the cells in this experiment
%     map_vec = iviv{folders};
%     %for all the desired cells
%     for cells = map_vec
%         
%         if mod(sub_count,25) == 1
%             h= figure;
%             sub_count = 1;
%         end
%         subplot(5,5,sub_count)
%         
%         %run the plotting function
%     %     sftf_plot(squeeze(sftf_all(cells,:,:,:)),sf_vec,tf_vec,dir_vec,h,squeeze(psth_all(cells,:,:,:)))
%         sftf_plot(squeeze(sftf_mat(folders).data(cells,:,:,:)),sf_vec,tf_vec,dir_vec,h)
%         
%         %update subplot counter
%         sub_count = sub_count + 1;
%     end
% end
%% Build a matrix with only the iviv cells

%allocate memory for the cells
sftf_all = zeros(sum(cellfun(@length,iviv)),length(sf_vec),length(tf_vec),length(dir_vec));
%initialize a counter
cell_count = 1;
%for all the folders
for folders = 1:folder_num
    %get the cells in this experiment
    map_vec = iviv{folders};
    %for all the desired cells
    for cells = map_vec
        %run the plotting function
        sftf_all(cell_count,:,:,:) = sftf_mat(folders).data(cells,:,:,:);
        %update the counter
        cell_count = cell_count + 1;
    end
end
%% Plot histograms of the values

close all
figure
%format the data to line all values in sf and tf coordinates only
% all_data = reshape(permute(sftf_all,[2 3 1 4]),3,3,[]);
%collapse direction dimension, then average across all cells
all_data = squeeze(nanmean(max(sftf_all,[],4),1));
imagesc(all_data)
set(gca,'XTick',[1 2 3],'XTickLabels',sf_vec,'YTick',[1 2 3],'YTickLabels',tf_vec,'YDir','normal')
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')


% plot a 2D histogram with the max from each cell
%alternatively, calculate the max across directions
prctile_data = squeeze(max(sftf_all,[],4));

%reformat matrix to take max
max_data = reshape(prctile_data,size(prctile_data,1),length(sf_vec)*length(tf_vec));
[~,sort_idx] = max(max_data,[],2);
%convert the idx into sf and tf
[thres_sf,thres_tf] = ind2sub([length(sf_vec),length(tf_vec)],sort_idx);
%plot the 2d histogram
figure
histogram2(thres_sf,thres_tf)
%label the axes
set(gca,'XTick',[1 2 3],'XTickLabels',sf_vec,'YTick',[1 2 3],'YTickLabels',tf_vec)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')

%plot a 2d histogram with the 80th percentile of the data

%then find the top 80th percentile and build a histogram
[sorted_data,sort_idx] = sort(prctile_data(:));
%get the percentile
prctile_val = prctile(sorted_data,80);
%get the values above threhold in the matrix, and the corresponding indexes
% thres_data = sorted_data(sorted_data>prctile_val);
thres_idx = sort_idx(sorted_data>prctile_val);
%convert the idx into sf and tf
[~,thres_sf,thres_tf] = ind2sub(size(prctile_data),thres_idx);

figure
histogram2(thres_sf,thres_tf)
%label the axes
set(gca,'XTick',[1 2 3],'XTickLabels',sf_vec,'YTick',[1 2 3],'YTickLabels',tf_vec)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')

%plot a histogram of the whole data distribution
figure
histogram(sftf_all(:))