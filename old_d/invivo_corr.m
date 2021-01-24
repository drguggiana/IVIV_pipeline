%% Clean up
clearvars
close all
%% Select the target files and load the ROIs

%define the load path
load_path = 'I:\Simon Weiler\IN VIVO_FINAL\spon activity\';

%pick the folders to use
folder_list = uipickfiles('FilterSpec',load_path);

%get the number of folders
exp_num = length(folder_list);
%allocate memory to store the ROI and neuropil data
ROI_exp = cell(exp_num,2);
%also for some metadata
meta_exp = cell(exp_num,1);
%for all the folders
for experiment = 1:exp_num
    
    %load the ROI data
    ROI_exp{experiment,1} = load(folder_list{experiment},'ROIs');
    ROI_exp{experiment,1} = ROI_exp{experiment,1}.ROIs;
    
    %load the neuropil data
    ROI_exp{experiment,2} = load(folder_list{experiment},'np');
    ROI_exp{experiment,2} = ROI_exp{experiment,2}.np;
    
    %load the relevant metadata
    meta_exp{experiment,1} = load(folder_list{experiment},'stimarray');
    meta_exp{experiment,1} = meta_exp{experiment,1}.stimarray;
    meta_exp{experiment,2} = load(folder_list{experiment},'info');
    meta_exp{experiment,2} = meta_exp{experiment,2}.info;
    
end
%% Preprocess each channel separately
close all
%allocate memory for the red and green channel data
f_cell = cell(exp_num,2);

%define the cell with the field names
field_names = {'activity','activity_r'};
%initialize a wait bar
h = waitbar(0,'Preprocessing progress');
%define the neuropil subtraction factor
r = 0.7;
%for all the experiments
for experiment = 1:exp_num
    %update the waitbar
    waitbar(experiment/exp_num,h)
    %for both channels (green first)
    for channel = 1:2
        %% Load the channel
        data = cat(1,ROI_exp{experiment,1}(:).(field_names{channel}))';
        npdata = cat(1,ROI_exp{experiment,2}(:).(field_names{channel}))';
        %get the dark frames
        dark_frames = ceil(max(meta_exp{experiment,2}.darkframes)/4);
        %trim both traces
        data = data(dark_frames+2:end,:);
        npdata = npdata(dark_frames+2:end,:);
        %% Lowpass the data at 0.8 Hz
        
        %the code below checks whether the frame rate was the same in the
        %previous iteration, so it doesn't recalculate a new filter every
        %time (makes it much faster)
        %if it's the first loop iteration
        if experiment==1&&channel==1
            %get the sampling rate (i.e. frame rate) Divided by four since
            %there are 4 planes
            sample_rate_old = meta_exp{experiment}.framerate/4;
            sample_rate = sample_rate_old;
        else %if not, 
            sample_rate = meta_exp{experiment}.framerate/4;
        end
        %if it's the first iteration or the frame rate changed
        if (sample_rate_old ~= sample_rate) || (experiment==1&&channel==1)
            %design the filter
            lpFilt = designfilt('lowpassfir','PassbandFrequency',0.8, ...
                     'StopbandFrequency',1,'StopbandAttenuation',65,'SampleRate',sample_rate);
        end
%         %visualize the filter
%         fvtool(lpFilt)

        %filter the data
        filt_data = filtfilt(lpFilt,data);
        filt_npdata = filtfilt(lpFilt,npdata);
%         %check the filtering
%         figure
%         plot((0:length(data)-1)./sample_rate,data(:,1))
% %         plot(data(:,1))
%         hold on
%         plot((0:length(data)-1)./sample_rate,filt_data(:,1))
%         plot(filt_data(:,1))
%         figure
%         imagesc(data)
        
        
        %% Subtract neuropil
        f_cell{experiment,channel} = filt_data-r.*filt_npdata+r.*median(filt_npdata,1);

    end
end
%close the waitbar
close(h)
%% Plot a target experiment for quality control
close all
%define the target experiment
tar_exp = 1;%:32;
%for all the target experiments
for exper = tar_exp
%     figure
%     imagesc(f_cell{exper,1})
    figure
    imagesc(f_cell{exper,2})
end
%% Calculate dRoR
close all

%define the percentile levels to explore
% prctile_levels = [1 10:10:100];
% prctile_levels = 80;
prctile_levels = 20;
%get the number of levels
level_num = length(prctile_levels);
%allocate memory to store the number of events detected for each percentile
events_mat = cell(level_num,1);

%for all the levels
for levels = 1:level_num
    
    fprintf(strcat('Current level: ',num2str(levels),'\r\n'))
    %allocate memory for the dRoR
    dRoR_cell = cell(exp_num,1);
    %initialize a wait bar
    h = waitbar(0,'dRoR calculation progress');
    %for all the experiments
    for experiment = 1:exp_num
        %update the waitbar
        waitbar(experiment/exp_num,h)
        %calculate R(t)
        R_t = f_cell{experiment,1}./f_cell{experiment,2};
        
        %detrend
        %get the sampling rate (i.e. frame rate)
        sample_rate = meta_exp{experiment}.framerate/4;
        %define the length in time of the window (in s)
        window_time = 14;
        %define the anonymous function to use in the moving window
        function_handle = @(x) prctile(x,8,1);
        
        %perform the detrending
        R_detrend = R_t - moving_window(R_t,round(sample_rate*window_time),function_handle,0);
        
        %calculate R0
        %define the length of the moving window (in s), from Winnubst et al.
        %2015
        r0_time = 20;
        %get the moving window function
        r0_function = @(x) prctile(x,prctile_levels(levels),1);
        %     r0_function = @(x) min(x,[],1);
        %     r0_function = @(x) mean(x,1);
        %calculate the baseline
%         R0 = median(moving_window(R_detrend,round(sample_rate*r0_time),r0_function));
        R0 = mean(moving_window(R_detrend,round(sample_rate*r0_time),r0_function,1));
        
        %calculate and save dRoR
        dRoR_cell{experiment} = (R_detrend-R0)./R0;
        %     %use R_detrend instead, the de-trended trace, since it is scaled better
        %     dRoR_cell{experiment} = R_detrend;
    end
    %close the waitbar
    close(h)
    %% OFF Plot a target experiment for quality control
%     close all
%     %define the target experiment
%     tar_exp = 2;
%     %define the target traces
%     tar_traces = 1:5;
%     %for all the target experiments
%     for exper = tar_exp
%     
%         figure
%         imagesc(((dRoR_cell{exper})'))
%     %     figure
%     %     histogram(dRoR_cell{tar_exp})
%         figure
%         plot((0:size(dRoR_cell{exper},1)-1)./sample_rate,dRoR_cell{exper}(:,tar_traces))
%     end
    %% Event detection
    
    %allocate memory for the results
    event_store = cell(exp_num,1);
    %for all the experiments
    for experiment = 1:exp_num
        %load the data
        data = dRoR_cell{experiment};
        %bin the data to 5Hz
        %define the binning factor
        bin_factor = 1.5;
        %allocate memory for the binned data
        bin_matrix = zeros(floor(size(data,1)/bin_factor),size(data,2));
        
        %get the binning map
        bin_map = discretize(1:size(data,1),1:bin_factor:size(data,1));
        
        %for all the bins
        for bins = 1:max(bin_map)
            %bin the data
            bin_matrix(bins,:) = mean(data(bin_map==bins,:),1);
        end
        %replace the original data matrix with the binned one for
        %processing
        data = bin_matrix;
        %zscore the data
%         data = zscore(data);
        %define the threshold
        %     thres = std(data,0,1);
%         thres = mean(data,1);
        thres = 2;
        %get the number of cells
        cell_num = size(data,2);
        %mark the places with supra-threshold signal 
        thres_map = [zeros(1,cell_num);zscore(diff(data,1,1))>thres];
        %get the sampling rate (i.e. frame rate)
        sample_rate = meta_exp{experiment}.framerate/4;
        %count the events
        events_percell = sum(thres_map,1)./(size(data,1)/(sample_rate/bin_factor));
%         %allocate memory to store the events per cell
%         events_percell = zeros(cell_num,1);
%         %get the event timings
%         [event_time_all,event_cell] = find(thres_map==1);
%         %for all the cells
%         for cells = 1:cell_num
%             %get the event times only for this cell
%             event_time = event_time_all(event_cell==cells);
%             %get the event ends
%             event_ends = event_time(diff(event_time)>1);
%             %and the event starts
%             event_starts = [event_time(1);event_time(find(diff(event_time(1:end))>1)+1)];
%             event_starts = event_starts(1:end-1);
% %             %get the event lenghts
% %             event_duration = event_ends-event_starts;
% %             %eliminate the events shorter than 2 frames and count
% %             events_percell(cells) = sum(event_duration>2);
%         end
        %store the results in the main cell
        event_store{experiment} = events_percell;
    end
%store the number of events
events_mat{levels} = event_store;
end
%% Plot results of the percentile finder

%if there are more than 1 level
if length(prctile_levels) > 1
    close all
    %allocate memory for the counts per animal
    counts_peranimal = zeros(exp_num,level_num);

    %for all the experiments
    for levels = 1:level_num
        %for all the experiments
        for experiment = 1:exp_num
            counts_peranimal(experiment,levels) = sum(vertcat(events_mat{levels}{experiment}));
        end
    end
    figure
    imagesc(counts_peranimal)
    xlabel('Percentile used')
    ylabel('Animal')
    title('Detected events as a function of animal and percentile')
    colorbar
    figure
    plot(counts_peranimal')
    xlabel('Percentile used')
    ylabel('Events detected')

    %% Define the percentile to use as the max

    %get the max for each animal
    [~,idx] = max(counts_peranimal,[],2);
    %set the event_store matrix as the mode of the max coordinate
    event_store = events_mat{mode(idx)};
end
%% Plot the results of the event detection
close all

figure
%for all the experiments
for experiment = 1:exp_num
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    histogram(event_store{experiment})
    xlabel('Firing rate')
end
%% Calculate the correlation metrics

close all
%allocate memory to store the correlation matrices
corr_mat = cell(exp_num,1);
figure
%for each cell
for experiment = 1:exp_num
%     %get the ROIs for this cell
%     ROIs = dRoR_cell{experiment};
    %load the traces into a matrix
    trace_mat = dRoR_cell{experiment};
    %calculate the correlation matrix for the activity traces
    corr_mat{experiment} = corr(trace_mat);
    %plot the matrix
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    imagesc(corr_mat{experiment})
end
%% Calculate population coupling

close all
%allocate memory to store the population couplings
popcop_cell = cell(exp_num,2);
figure
%for each cell
for experiment = 1:exp_num
%     %get the ROIs for this cell
%     ROIs = ROI_exp{experiment};
    %load the traces into a matrix
    trace_mat = dRoR_cell{experiment};
    %get the number of traces
    trace_num = size(trace_mat,2);
%     %calculate the average signal
%     ave_signal = mean(trace_mat,2);
    
    %allocate memory for the pop couplings
    pop_mat = zeros(trace_num,1);
    %for all the traces
    for traces = 1:trace_num
        %get the target trace
        tar_trace = trace_mat(:,traces);
        %get a vector pointing to the other traces
        point_vec = ones(trace_num,1)==1;
%         point_vec(traces) = 0;
        %get the signal from the rest of the cells
        rest_signal = mean(trace_mat(:,point_vec),2);
        %calculate the correlation
        pop_mat(traces) = corr(tar_trace,rest_signal);
    end
%     %calculate the correlation matrix for the activity traces
%     popcop_cell{experiment} = corr(ave_signal,trace_mat);
    %store the matrix in the main cell
    popcop_cell{experiment,1} = zscore(pop_mat);
    %also store the number of traces going into the calculation
    popcop_cell{experiment,2} = trace_num;
    %plot the matrix
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    histogram(popcop_cell{experiment,1})
    title(num2str(trace_num))
end

%plot the overall population coupling histogram
figure
histogram(vertcat(popcop_cell{:,1}))
%% Calculate the "noise" correlation

close all
%allocate memory to store the correlation matrices
noise_mat = cell(exp_num,1);
figure
%for each cell
for experiment = 1:exp_num
%     %get the ROIs for this cell
%     ROIs = ROI_exp{experiment};
    %load the traces into a matrix
    trace_mat = dRoR_cell{experiment};
    %subtract the mean from each trace
    trace_mat = trace_mat-mean(trace_mat,2);
    %calculate the correlation matrix for the activity traces
    noise_mat{experiment} = corr(trace_mat);
    %plot the matrix
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    imagesc(noise_mat{experiment})
end
%% Save the calculated metrics

%define the save path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Spont_activity_out';

%assemble the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_spontAct.mat');

%save the file
% save(fullfile(save_path,save_name),'popcop_cell','event_store','noise_mat',...
%     'corr_mat','prctile_levels','dRoR_cell')