function sftf_plot(data_in,varargin)

%this function will take a 3D matrix and plot it. The rows will be plotted as
%spatial frequency, the columns as temporal frequency and the preferred
%direction as the direction of the arrow patch

%input list
%1) data matrix, as a 3D matrix with SF, TF and direction as dimensions
%2) vector with the spatial frequencies used
%3) vector with the temporal frequencies used
%4) vector with the direction angles used
%5) handle of the target figure
%6) OPTIONAL psth averages to plot

%get the dimensions of the matrix
sf_dim = size(data_in,1);
tf_dim = size(data_in,2);
% dir_dim = size(data_in,3);
%get the target figure handle
h = varargin{4};

%get the labels for sf, tf and the vector for dir
sf_label = varargin{1};
tf_label = varargin{2};
dir_vec = varargin{3};
%define the color code and depth
c_map = parula(256);
% dir_vec = [0 45 90 135 180 225 270 315]
%% Calculate the max

%allocate memory to store the calcium max integral and the pref direction
plot_matrix = zeros(sf_dim,tf_dim,2);
%for all the sf
for sf = 1:sf_dim
    %for all the tf
    for tf = 1:tf_dim
        %get the max value and direction for this combination of sftf
        [plot_matrix(sf,tf,1),plot_matrix(sf,tf,2)] = max(data_in(sf,tf,:)); 
    end
end
%% Plot the matrix

%target the provided figure
figure(h)

% subplot(1,30,1:28)
%define the patch coordinates
% x_data = [0 1 2 1.5 1.5 0.5 0.5]';
% y_data = [1 2 1 1 0 0 1]';
x_data = [1 0 1 1 2 2 1];
y_data = [0 1 2 1.5 1.5 0.5 0.5];
p_center = [1,1,0];
%define the arrow distance factor
dist_f = 3;
%for all the sf
for sf = 1:sf_dim
    %for all the tf
    for tf = 1:tf_dim
%         subplot(sf_dim,tf_dim,sub_count)
        %if the value is a NaN, skip the iteration
        if isnan(plot_matrix(sf,tf,1))
            continue
        end
%         %normalize to cell
%         plot_matrix(:,:,1) = normr_2(plot_matrix(:,:,1)).*255;
        %create the patch
        p = patch(x_data+sf*dist_f,y_data+tf*dist_f,c_map(floor(plot_matrix(sf,tf,1))+1,:));
        rotate(p,[0 0 1],-dir_vec(plot_matrix(sf,tf,2)),p_center + [sf*dist_f,tf*dist_f,0])
        
        
    end
end

%set axes properties
axis square
box on
set(gca,'XTick',[],'YTick',[],'TickLength',[0 0])
set(gca, 'color', 'none');
set(gca,'XLim',[2 12],'YLim',[2 12])
%label the axes
set(gca,'XTick',[4 7 10],'XTickLabels',sf_label,'YTick',[4 7 10],'YTickLabels',tf_label)
xlabel('Spatial Frequency (cyc/deg)')
ylabel('Temporal Frequency (cyc/s)')

%plot a colorbar
colorbar

% %plot the colorbar
% subplot(1,30,29:30)
% 
% imagesc((0:255)')
% xlabel('A.U.')
% set(gca,'XTick',[],'YTick',[1 255],'TickLength',[0 0],'YTickLabels',[0 1],...
%     'ydir','normal')
%% Optional psth plotting

%check the number of inputs. if enough, plot the psth too
if nargin > 5
    
    %load the psth info
    psth_traces = varargin{5};
    figure

    %define the amplitude factor
%     amp_rat = 2;
    amp_rat = 0.5;
    %define the subsampling factor
%     sub_rat = 10;
    sub_rat = 100;
    %define the separation factor
    sep_factor = (size(psth_traces,3) + 1000)/sub_rat;
    %for all the traces
    for x = 1:size(plot_matrix,1)
        for y = 1:size(plot_matrix,2)
            
            %select the trace color depending on the response type
            trace_resp = isnan(plot_matrix(x,y,1));
            switch trace_resp
                case 0
                    trace_color = 'b';
                case 1
                    trace_color = 'r';
                case 2
                    trace_color = 'm';
                otherwise
                    trace_color = 'k';
            end

            %get the x vector, correcting for the array position
%             x_vec = (1:length(squeeze(psth_traces{plot_matrix(sf,tf,2),x,y}))) + sep_factor*x;
            %and the y vector
%             y_vec = squeeze(psth_traces{plot_matrix(sf,tf,2),x,y})./amp_rat + sep_factor*y;

                        %calculate the tuning curve of the sftf combo and plot
                        tc = normr_2(-tc_calc(psth_traces(:,x,y)));
                        x_vec = (1:length(tc)) + sep_factor*x;
                        y_vec = (tc./amp_rat) + sep_factor*y;

%             %for all directions
%             for dirs = 1:size(psth_traces,1)
%                 y_vec = squeeze(psth_traces{dirs,x,y})./amp_rat + sep_factor*y;
% 
                %plot the result
                plot(x_vec,y_vec,trace_color)
                hold('on')
%             end
            %plot a 0 line
            plot(x_vec([1 end]),[0 0]+ sep_factor*y,'k--')
%             %plot a line at 7 ms
%             plot(x_vec([7 7]),[max(y_vec) min(y_vec)],'g--')
            
%             %plot the cluster of origin
%             text(x_vec(1),-y_vec(1),num2str(clu_idx(trace_id(:,1)==curr_ind_clu &trace_id(:,2)==target_map)))
            
            
        end
    end
    %configure the plot
    set(gca,'YLim',[8 35])
end
function tc = tc_calc(psth_in)

%assemble the matrix for svd
svd_mat = cat(1,psth_in{:});

%turn NaNs to 0
svd_mat(isnan(svd_mat)) = 0;

%run the svd
[U,~,~] = svd(normr_2(svd_mat));

%output the tuning curve
tc = U(:,1);