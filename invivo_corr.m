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
% %         plot(filt_data(:,1))
%         figure
%         imagesc(data)
        
        
        %% Subtract neuropil
        f_cell{experiment,channel} = filt_data-r.*filt_npdata+r.*median(filt_npdata,1);

    end
end

%Plot a target experiment for quality control
%define the target experiment
tar_exp = 19;
figure
imagesc(f_cell{tar_exp,1})
figure
imagesc(f_cell{tar_exp,2})
%close the waitbar
close(h)
%% Calculate dRoR
close all
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
    R_detrend = R_t - moving_window(R_t,round(sample_rate*window_time),function_handle);
    
%     %calculate R0
%     %define the length of the moving window (in s)
%     r0_time = 20;
%     %get the moving window function
%     r0_function = @(x) prctile(x,8,1);
%     r0_function = @(x) min(x,[],1);
%     r0_function = @(x) mean(x,1);
%     %calculate the baseline
%     R0 = median(moving_window(R_detrend,round(sample_rate*r0_time),r0_function));
    
%     %calculate and save dRoR
%     dRoR_cell{experiment} = (R_detrend-R0)./R0;
    %use R_detrend instead, the de-trended trace, since it is scaled better
    dRoR_cell{experiment} = R_detrend;
end

%Plot a target experiment for quality control
%define the target experiment
tar_exp = 6;
figure
imagesc(((dRoR_cell{tar_exp})'))
figure
histogram(dRoR_cell{tar_exp})
figure
plot(dRoR_cell{tar_exp}(:,1:5))
%close the waitbar
close(h)
%% Event detection

%allocate memory for the results
event_store = cell(exp_num,1);
%for all the experiments
for experiment = 1:exp_num
    %load the data
    data = dRoR_cell{experiment};
    %define the threshold
%     thres = std(data,0,1);
    thres = mean(data,1);
    %mark the places with supra-threshold signal
    thres_map = data>2.*thres;
    %get the number of cells
    cell_num = size(thres_map,2);
    %allocate memory to store the events per cell
    events_percell = zeros(cell_num,1);
    %get the event timings
    [event_time_all,event_cell] = find(thres_map==1);
    %for all the cells
    for cells = 1:cell_num
        %get the event times only for this cell
        event_time = event_time_all(event_cell==cells);
        %get the event ends
        event_ends = event_time(diff(event_time)>1);
        %and the event starts
        event_starts = [event_time(1);event_time(find(diff(event_time(1:end))>1)+1)];
        event_starts = event_starts(1:end-1);
        %get the event lenghts
        event_duration = event_ends-event_starts;
        %eliminate the events shorter than 2 frames and count
        events_percell(cells) = sum(event_duration>2);
    end
    %store the results in the main cell
    event_store{experiment} = events_percell;
end
%plot the results
close all

figure
%for all the experiments
for experiment = 1:exp_num
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    histogram(event_store{experiment})
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
popcop_cell = cell(exp_num,1);
figure
%for each cell
for experiment = 1:exp_num
%     %get the ROIs for this cell
%     ROIs = ROI_exp{experiment};
    %load the traces into a matrix
    trace_mat = dRoR_cell{experiment};
    
    %calculate the average signal
    ave_signal = mean(trace_mat,2);

    %calculate the correlation matrix for the activity traces
    popcop_cell{experiment} = corr(ave_signal,trace_mat);
    %plot the matrix
    subplot(ceil(sqrt(exp_num)),round(sqrt(exp_num)),experiment)
    histogram(popcop_cell{experiment})
end

%plot the overall population coupling histogram
figure
histogram(horzcat(popcop_cell{:}))
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