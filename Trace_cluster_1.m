%% Clean up
clearvars
close all
addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load paths

%pick the folders to use
folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\INPUT MAPS_final\Setup2');


%define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';
%get the number of folders
folder_num = length(folder_list);
%allocate memory to store the subfolders that qualify
folder_cell = cell(folder_num,1);
%for all the folders
for folders = 1:folder_num

    %get the subfolders in the first level of this folder
    subfolder_list = dir(folder_list{folders});

    %also, get rid of the non-folders
    subfolder_list = subfolder_list(vertcat(subfolder_list(:).isdir)==1);
    %leave only the folders that have SW in them
    subfolder_list = subfolder_list(contains({subfolder_list(:).name},'SW'));
    %finally, build the full paths for the folders, including the "images"
    %folder and the "map01"
    path_list = cell(size(subfolder_list,1),1);
    %for all the subfolders
    for subs = 1:size(subfolder_list,1)

        %load all the map paths in a given cell
        
        %get the list of paths
        xsg_paths = dir(strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\*map*'));
        %and create a cell to store the full paths
        xsg_cell = cell(length(xsg_paths),1);
        %also a vector to keep track of the empty folders
        elim_vec = ones(size(xsg_cell,1),1);
        %for each path
        for paths = 1:length(xsg_paths)
            %load the files in the directory
            xsg_files = dir(strcat(xsg_paths(paths).folder,'\',xsg_paths(paths).name,'\*.xsg'));
            %if there is no file, skip the entry and update the vector
            %accordingly
            if isempty(xsg_files)
                elim_vec(paths) = 0;
                continue
            end
            %append the name of the file to the path
            xsg_cell{paths} = strcat(xsg_files(1).folder,'\',xsg_files(1).name);
        end
        %store the cell with paths in the larger storage cell
        path_list{subs} = xsg_cell(elim_vec == 1);
    end
%     %update the folder cell only with the cells that had an image
%     folder_cell{folders} = path_list(elim_vec==1,:);
    %update the folder cell
    folder_cell{folders} = vertcat(path_list{:});
end

%concatenate the entire list of paths
folder_all = vertcat(folder_cell{:});
%get the number of maps to do
map_num = size(folder_all,1);
%% Get an array with the list of cell names for each map

%extract the cell names from the folder all array
%allocate memory to save the cell names
cell_names = cell(map_num,1);
%for all the folders
for maps = 1:map_num
    temp_name = strsplit(folder_all{maps},'\');
    cell_names{maps} = strcat(temp_name{end-2}(1:2),temp_name{end-3}(1:6),temp_name{end-2}(3:end));
end
%% Load, trim and reorder the data for the PCA 

%load a random 10% of the data and run PCA on it
%get the random indexes of the maps to load
% rand_ind = randperm(map_num,round(map_num/10));
% rand_ind = randperm(map_num,68);
%allocate memory to store the maps
temp_xsg = load(folder_all{1},'-mat');
%get the size of the maps
map_size = length(temp_xsg.data.ephys.trace_1);
%also the number of positions
num_positions = temp_xsg.header.mapper.mapper.positionNumber;
%define the interval to grab. Each recording is a second and the desired
%range is from 100 to 300 ms
trace_range = 1000:3000;
%define the time interval for the background
trace_background = 1:999;

%allocate a matrix for them
% map_matrix = zeros(length(trace_range),length(rand_ind));
map_matrix = zeros(length(trace_range)*num_positions,map_num);
%allocate memory for the standard deviation of the background for each
%trace
background_std = zeros(num_positions,map_num);

%set up a vector for the skipped files
elim_vec = ones(size(map_matrix,2),1);

%create a waitbar
w_bar = waitbar(0,'Loading maps');
%load the random maps into the matrix
% for maps = 1:length(rand_ind)
for maps = 1:map_num
    
    %show progress
    waitbar(maps/map_num,w_bar)
    %load the xsg file
%     temp_xsg = load(folder_all{rand_ind(maps)},'-mat');
    temp_xsg = load(folder_all{maps},'-mat');
    %check the number of grid points in the cell
    temp_map = temp_xsg.header.mapper.mapper.mapPatternArray;
    %if the trace has less than 16x16
    if size(temp_map,1) < 16 || size(temp_map,2) < 16 ||...
            isempty(temp_xsg.data.ephys)
        %mark the position and skip the trace
        elim_vec(maps) = 0;
        continue
    end
    %load the time trace
    temp_trace = temp_xsg.data.ephys.trace_1;
    %reshape to trim the relevant part of the data only
    temp_trace = reshape(temp_trace,[],num_positions);
   
    %get the background activity
    background_act = mean(temp_trace(trace_background,:),1);
    %get the std deviation of the background
    background_std(:,maps) = std(temp_trace(trace_background,:),0,1);
    %trim the trace
    temp_trace = temp_trace(trace_range,:);
    %subtract background activity
    temp_trace = temp_trace - background_act;
    %also load the map order to have all the maps ordered in their location
    %instead of in time of stimulation
    %linearize the map and apply to the data
    temp_trace = temp_trace(:,temp_map(:));
    %reshape again and store
    map_matrix(:,maps) = reshape(temp_trace,length(trace_range)*num_positions,[]);
end

close(w_bar)

%get rid of the empty spaces with the skipped traces
map_matrix = map_matrix(:,elim_vec==1);
%also within the std matrix
background_std = background_std(:,elim_vec==1);
% map_matrix2 = map_matrix;
%and also modify the folder vector (so I can refer back to the particular
%maps)
folder_all = folder_all(elim_vec==1);
%update the cell name array
cell_names = cell_names(elim_vec==1);

%get the unique cell names
[uni_cells,~,cell_id] = unique(cell_names);
%% Split the maps into excitatory and inhibitory

%calculate the average value for each map, since area under the curve
%should report whether it is exc or inh
ave_value = mean(map_matrix,1);

figure
histogram(ave_value)
%calculate the polarity of each map to propagate later (0 is exc)
map_polarity = ave_value>0;
%% Reshape the matrix, so that all traces are in a single column

map_matrix2 = reshape(map_matrix,length(trace_range),[]);
%linearize the std matrix also
background_std = background_std(:);
%and create a map from traces to folders (including polarity)
%allocate memory for the map
trace2folder = zeros(num_positions,size(map_matrix,2),5);
%for all the folders
for folders = 1:size(folder_all,1)
    trace2folder(:,folders,1) = 1:num_positions;
    trace2folder(:,folders,2) = folders;
    trace2folder(:,folders,3) = map_polarity(folders);
    trace2folder(:,folders,4) = cell_id(folders);
end
%also reshape this map to be able to refer to the original map
trace2folder = reshape(trace2folder,[],5);
%% Filter out traces that are too flat (using 3 std criterion)
close all
%calculate the std of each trace
std_matrix = std(map_matrix2,0,1);
% %define the percentile threshold
% prctile_thres = 80;
% prctile_cutoff = prctile(std_matrix,prctile_thres);
%plot the distribution of the std
figure
histogram(std_matrix)
% %exclude traces on the lowest 10th percentile
% map_matrix2 = map_matrix2(:,std_matrix>prctile_cutoff);
%exclude traces with less than 3 stds over background
map_matrix2 = map_matrix2(:,std_matrix'>3.*background_std);
%also save a map of the selected traces
prctile_map = std_matrix'>3.*background_std;
%generate a map of the origin of the traces left after the percentile cut
postperc_map = trace2folder(prctile_map,:);
%mark the places on the trace2folder matrix
trace2folder(prctile_map==0,5) = NaN;
%% Bin the data

%temporally bin the data by a defined factor
%define the factor
bin_factor = 10;

%bin by the factor

%allocate memory for the binned data
bin_matrix = zeros(floor(size(map_matrix2,1)/bin_factor),size(map_matrix2,2));

%get the binning map
bin_map = discretize(1:size(map_matrix2,1),1:bin_factor:size(map_matrix2,1));

%for all the bins
for bins = 1:max(bin_map)
    %bin the data
    bin_matrix(bins,:) = mean(map_matrix2(bin_map==bins,:),1);
end
%% OFF Split the data in inhibitory and excitatory traces

% close all
% %allocate memory for the trace matrices
% separate_maps = cell(2,1);
% 
% %calculate the average value for each trace, since area under the curve
% %should report whether it is exc or inh
% ave_value = mean(bin_matrix,1);
% 
% figure
% histogram(ave_value)
% %split the traces (1 is exc, 2 is inh)
% separate_maps{1} = bin_matrix(:,ave_value<0);
% separate_maps{2} = bin_matrix(:,ave_value>0);
%% Analyze exc and inh separately

%allocate memory for the results
EI_cell = cell(14,2);
%create a waitbar
h = waitbar(0,'Clustering progress');
%for exc and inh
for exc_inh = 1:2
    %% Run the PCA
    %run a pca on the data
    [coeff,score,latent] = pca((bin_matrix(:,postperc_map(:,3)==exc_inh-1)'));

    close all
    %plot the relevant outputs
    figure
    plot(latent)
    hold('on')
    plot(cumsum(latent)/sum(latent),'k--')

    figure
    imagesc(coeff)
    %% Skip the first PC and use the other ones to cluster

    %define the PCs to use
    pc_vec = 1:21;
    %get the score matrix with the PCs to use
    clu_mat = (score(:,pc_vec));

    %cluster the matrix using GMMs

    %define a vector with cluster numbers
    % clu_vec = [5 10 15 20 50 100];
    clu_vec = [10 15 20 25 30 35 40 50 100 200];
    % clu_vec = 50;
    %get the number of cluster runs
    clu_size = length(clu_vec);
    %create a cell array to save the clustering results
    clu_cell = cell(clu_size,3);
    %define the statistical settings
    opts = statset('MaxIter',1000);
    %define the number of replicates
    rep_num = 10;
    %for all the cluster numbers
    for clu = 1:clu_size
        waitbar((clu + clu_size*(exc_inh-1))/(clu_size*2),h)
%         %show the current cluster number
%         fprintf(strcat('Current cluster number:',num2str(clu_vec(clu)),'\r\n'))
        %create a GMM for the data
        clu_cell{clu,1} = fitgmdist(clu_mat,clu_vec(clu),'RegularizationValue',0.0001,...
            'CovarianceType','diagonal','Replicates',rep_num,'Options',opts);

        %cluster the data accordingly
        clu_cell{clu,2} = cluster(clu_cell{clu,1},clu_mat);

        %and get the BIC value
        clu_cell{clu,3} = clu_cell{clu,1}.BIC;
    end

    %plot the BIC
    figure
    plot(clu_vec,vertcat(clu_cell{:,3}));

    %get the BIC minimum coordinate
    [~,bic_min] = min(vertcat(clu_cell{:,3}));
    %and the associated number of morpho_clusters
    clunum = clu_vec(bic_min);
    %get the indexes from the best model
    clusters = clu_cell{bic_min,2};
    %% Plot the clustering results
    close all

    %define which matrix to use for the traces
    plot_matrix = map_matrix2(:,postperc_map(:,3)==exc_inh-1);
    %also get the map to original traces only for exc/inh
    temp_postperc_map = postperc_map(postperc_map(:,3)==exc_inh-1,:);

    %allocate memory to store the cluster averages
    clu_ave = zeros(clunum,size(clu_mat,2));
    %save the number of members too
    clu_mem = zeros(clunum,1);
    %and an average of the traces going into each cluster
    trace_ave = zeros(clunum,size(plot_matrix,1));
    trace_std = zeros(clunum,size(plot_matrix,1));

    %for all the morpho_clusters
    for clu = 1:clunum
        %average the cells in question
        clu_ave(clu,:) = squeeze(mean(clu_mat(clusters==clu,:),1));
        %count the number of members
        clu_mem(clu) = sum(clusters==clu);
        %and average the actual traces
        trace_ave(clu,:) = mean(plot_matrix(:,clusters==clu),2);
        trace_std(clu,:) = std(plot_matrix(:,clusters==clu),0,2)./sqrt(clu_mem(clu));
    end

    figure
    imagesc(clu_ave)

    %also plot a cluster average of the actual traces in each cluster
    figure
    imagesc(trace_ave)

    %generate a colormap for the cluster averages
    c_map = jet(clunum);
    %define the length of the trace to plot
    plot_length = 1:2001;
    figure
    %for all the clusters
    for clu = 1:clunum
%         figure
    %     plot(1:size(plot_matrix,1),trace_ave(clu,:)'+clu)
        %define the x vector
        x_vec = (0:size(plot_length,2)-1)./10;
        subplot(round(sqrt(clunum)),ceil(sqrt(clunum)),clu)
        shadedErrorBar(x_vec,trace_ave(clu,plot_length)'...
            ,trace_std(clu,plot_length)','lineprops',{'Color',c_map(clu,:)})
        hold('on')
        %plot the zero line
        plot(x_vec',zeros(length(x_vec),1),'k--')
        plot([7 7],get(gca,'YLim'),'k--')
        ylabel('Current (pA)')
        xlabel(strcat('Time (ms) #:',num2str(clu_mem(clu))))

    end
    %store the results in the main cell
    EI_cell{1,exc_inh} = clu_cell;
    EI_cell{2,exc_inh} = clu_mat;
    EI_cell{3,exc_inh} = trace_ave;
    EI_cell{4,exc_inh} = trace_std;
    EI_cell{5,exc_inh} = x_vec;
    EI_cell{6,exc_inh} = clu_vec;
    EI_cell{7,exc_inh} = opts;
    EI_cell{8,exc_inh} = clunum;
    EI_cell{9,exc_inh} = pc_vec;
    EI_cell{10,exc_inh} = coeff;
    EI_cell{11,exc_inh} = score;
    EI_cell{12,exc_inh} = latent;
    EI_cell{13,exc_inh} = clusters;
    EI_cell{14,exc_inh} = clu_mem;    
end
close(h)
%% Look up table for the clusters

%define the code to use
%0: throw away
%1: synaptic
%2: direct AMPA
%3: direct AMPA + NMDA
%4: mixed 

%allocate memory for the table
EI_lut = cell(2,1);

%zscored
% EI_lut{1} = [4,2,3,2,2,2,0,4,4,3,...
%     2,2,0,0,4,1,4,3,3,2,...
%     4,2,4,2,2,4,4,2,3,4,...
%     1,2,3,2,3,3,3,4,2,2,...
%     1,2,3,1,1,3,4,4,1,4];

%non zscored
EI_lut{1} = [2,2,2,4,2,2,2,1,4,1,...
    2,1,3,3,3,1,2,3,3,2,...
    1,1,3,2,3,4,2,2,4,4,...
    3,3,2,1,1,3,2,1,4,4,...
    1,1,2,1,2,2,4,4,2,0];

%0: throw away
%1: synaptic
%2: direct AMPA
%3: delayed
%4: mixed direct and synaptic

%zscored
% EI_lut{2} = [2,1,1,1,4,1,1,1,0,4,...
%     1,1,4,2,1,1,1,1,4,2,...
%     1,0,4,2,1,1,1,1,2,4,...
%     4,1,1,4,1,2,4,4,1,4,...
%     1,4,1,0,4,1,1,1,4,1];

%non zscored
EI_lut{2} = [4,1,1,4,3,4,1,4,4,1,...
    1,4,1,2,2,4,4,3,4,4,...
    1,4,4,1,4,2,1,4,1,1,...
    1,1,4,4,1,1,4,4,4,4,...
    4,3,1,1,4,1,4,3,1,1];
%% Convert the cluster numbers into response type numbers

%for excitation and inhibition
for polarity = 1:2
    %load the cluster vector
    clusters = EI_cell{13,polarity};
    
    %change the values to their corresponding response types
    resp_types = EI_lut{polarity}(clusters);
    
    %update the trace2folder matrix with the response types
    trace2folder(prctile_map&trace2folder(:,3)==polarity-1,5) = resp_types;
end
%% Save the cluster data

% %define the save path
% save_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';

%define the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_rawClusters.mat');

%save the clustering matrix and the cluster indexes
save(fullfile(out_path,save_name),'folder_all','rep_num','prctile_map',...
    'trace_background','trace_range','EI_cell','EI_lut','uni_cells',...
    'trace2folder','background_std')
error('Done')
%% Determine the number of response types per map

close all

%get the map to original traces only for exc/inh
temp_postperc_map = trace2folder(trace2folder(:,3)==1,:);

figure
%show the maps that are present at this stage at all
%create a vector to show the presence of each map
map_vector = zeros(size(map_matrix,2),1);
%fill in the values present in the filtered data set
map_vector(unique(temp_postperc_map(:,2))) = 1;
imagesc(map_vector)
colormap(jet)
ylabel('Map #')
title('Presence of a given map after filtering by std percentile')

figure
%allocate a matrix with maps and clusters
map2clu = zeros(size(map_matrix,2),clunum);
%for all the clusters
for clu = 1:clunum
    map2clu(unique(temp_postperc_map(clusters==clu,2)),clu) = 1;
end

imagesc(map2clu)
colormap(jet)
ylabel('Map #')
xlabel('Clusters')
title('Presence of a particular response type within a map')
%% Pull the full 256 locations for a given map (and mark the ones for a given cluster)
close all

%define the cluster to focus on
resp_type = [1 4];
%set whether to look at exc or inh
exc_inh = 0;
%get the IDs of the traces corresponding to this cluster
% trace_id = temp_postperc_map(clusters==clu_focus,:);
trace_id = trace2folder(any(trace2folder(:,5)==resp_type,2)&trace2folder(:,3)==exc_inh,:);
%get the map ids for all traces in this cluster
map_ids = unique(trace_id(:,2));

%for all the maps
for target_map = map_ids(randperm(length(map_ids),20))'
    %define the map of interest
%     target_map = 60;
    figure
    %get the map from the main matrix
    map_traces = permute(reshape(map_matrix(:,target_map),[],sqrt(num_positions),sqrt(num_positions)),[3 2 1]);

    %define the amplitude factor
    amp_rat = 2;
    %define the subsampling factor
    sub_rat = 10;
    %define the separation factor
    sep_factor = (size(map_traces,3) + 1000)/sub_rat;
    %for all the traces
    for x = 1:size(map_traces,1)
        for y = 1:size(map_traces,2)
            %get the index corresponding to this position in single index
            curr_ind = sub2ind([size(map_traces,1),size(map_traces,2)],x,y);
            curr_ind_clu = sub2ind([size(map_traces,1),size(map_traces,2)],y,x);
            %if the trace belongs to the target color
            if any(trace_id(:,1)==curr_ind_clu & trace_id(:,2)==target_map)
                %make it red
                trace_color = 'r';
            else
                %otherwise make it blue
                trace_color = 'b';
            end
            %get the x vector, correcting for the array position
            x_vec = (1:length(squeeze(map_traces(x,y,1:sub_rat:end)))) + sep_factor*x;
            %and the y vector
            y_vec = squeeze(map_traces(x,y,1:sub_rat:end))./-amp_rat + sep_factor*y;
            %plot the result
            plot(x_vec,-y_vec,trace_color)
            hold('on')
            %plot a 0 line
            plot(x_vec([1 end]),[0 0]- sep_factor*y,'k--')
        end
    end
end
%% Use ICA to process the traces with a direct and a synaptic response

close all

%define the target trace
tar_clu = 2;
%define the number of components
num_components = 3;
%show the target trace
plot((0:size(plot_matrix,1)-1)./10,trace_ave(tar_clu,:))

%get the traces for this cluster
clu_traces = plot_matrix(:,clusters==tar_clu);

%use rica
mdl = rica(clu_traces,num_components,'Standardize',true);

% tranf_comp = transform(mdl,trace_ave(tar_trace,:));

transf_weights = mdl.TransformWeights;

figure
% plot(clu_traces(1,:)')
% hold('on')
%for all the components
for comps = 1:num_components
    plot(transf_weights(:,comps))
    
    hold('on')
%     plot(transf_weights(:,comps).*clu_traces(1,:)','--')
end
%% Back transform the traces into the demixed data
close all

figure
%back transform the data
transf_data = transform(mdl,clu_traces);
%plot it
plot((0:size(plot_matrix,1)-1)./10,transf_data)
%% Use Fourier spectral decomposition to split the traces

close all

%specify the parameters
Fs = 10000; %10kHz sampling
T = 1/Fs;
L = 2001;
t = (0:L-1)*T;
%define the target trace
tar_clu = 9;

figure
%show the target trace
plot((0:size(plot_matrix,1)-1)./10,trace_ave(tar_clu,:))

% %use fft to get the spectrum of the average trace
% fft_trace = fft(trace_ave(tar_clu,:));

% %compute spectrum
% %define the amplitudes from the raw fft
% %first take the real part of the output, scaled by the duration
% P2 = abs(fft_trace/L);
% %now take only one half of it (don't need it mirrored)
% P1 = P2(1:round(L/2)+1);
% %and rescale the middle (no idea why here :p)
% P1(2:end-1) = 2*P1(2:end-1);
% %define the frequency vector
% f = Fs*(0:(round(L/2)))/L;
% figure
% %plot the spectrum
% plot(f,P1)
%% Compute the spectrum for all the traces

close all

%get the traces for the cluster defined above
clu_traces = plot_matrix(:,clusters==tar_clu)';
%get the number of traces
clu_tracenum = size(clu_traces,1);
% %allocate memory to store the spectra
% spec_mat = zeros(clu_tracenum,round(size(clu_traces,2)/2));

spec_mat = fft(clu_traces,[],2);

%first take the real part of the output, scaled by the duration
P2 = abs(spec_mat/L);
%now take only one half of it (don't need it mirrored)
P1 = P2(:,1:round(L/2)+1);
%and rescale the middle (no idea why here :p)
P1(:,2:end-1) = 2*P1(:,2:end-1);
%define the frequency vector
f = Fs*(0:(round(L/2)))/L;

% figure
% imagesc(P1)

%plot all the spectra
figure
plot(f,P1)
set(gca,'XScale','log')

%plot the average of all the spectra
figure
shadedErrorBar(f,mean(P1),std(P1))
set(gca,'XScale','log')
%% Try GMMs for splitting the traces

close all

%define the target trace
tar_clu = 2;

%define the number of components
num_components = 3;

%fit the model
GMM = fitgmdist();
