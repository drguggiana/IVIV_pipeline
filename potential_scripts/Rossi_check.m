% Test the slope hypothesis

%% Clean up

clearvars
close all
matlabrc
Paths
%% Load the structure

% load the main structure
main_path = stage5_fullIVIV_path;
str = load(find_newer_file(main_path));
str = str.str;
%% Get the conversion between visual space and cortical space

% conversion factor from Garrett et al from mm to deg
conversion_factor = 36; %(um/deg)

% distance for the input map
input_distance = 69*15/2;

% angle in visual space
visual_deg = input_distance/conversion_factor;

% define the distance from mouse (cm)
mouse_distance = 100;

% distance in visual space (cm)
visual_distance = 2*pi*mouse_distance*(visual_deg/360);

% define the axis labels
x_axis = linspace(0,visual_deg,8);

% from Rossi
min_delta = 0.8; %(deg)
max_delta = 2.4; %(deg)

time_peak = 10; %(ms)

% speed of an object moving that distance (assume 2 divisions in inputs map from base to peak)
speed_object = 2*69/(time_peak*0.001);
%% Parameter plot

% define the vector of distance to object(cm)
object_vector = 1:100;

% define the vector of distance travelled in the slice (um)
slice_vector = 1:200;

% asume the time (ms)
target_time = 10;

% allocate memory for the results
speed_map = zeros(length(object_vector),length(slice_vector));
% start counters for both variables
d_count = 1;
% for all the object distances
for distance = object_vector
    t_count = 1;

    % for all the travel distances
    for travel = slice_vector
        % get the angle
        angle_travel = travel/conversion_factor;
        % get the distance travelled in visual space
        visual_distance = 2*pi*distance*(angle_travel/360);
        % get the speed and save
        speed_map(d_count,t_count) = visual_distance/(target_time*1e-3);
        % update the counter
        t_count = t_count + 1;
    end
    % update the counter
    d_count = d_count + 1;
end

% Plot
close all
figure

% define the density of labels
dens = 10;
imagesc((speed_map'))

set(gca,'TickLength',[0 0],'XTick',1:dens:length(object_vector),'XTickLabel',object_vector(1:dens:end),'XTickLabelRotation',0)
set(gca,'YTick',1:dens:length(slice_vector),'YTickLabel',slice_vector(1:dens:end))
xlabel('Object distance (cm)')
ylabel('Slice distance (um)')

% superimpose the Rossi limits
hold on
min_dist = min_delta*conversion_factor;
plot(get(gca,'XLim'),[min_dist min_dist].*max(slice_vector)./length(slice_vector),'c')

max_dist = max_delta*conversion_factor;
plot(get(gca,'XLim'),[max_dist max_dist].*max(slice_vector)./length(slice_vector),'r')

set(gcf,'Color','w')
h = colorbar;
ylabel(h,'Speed (cm/s)')
colormap(viridis)
%% Plot with only the angle dependency

% define the vector of distance travelled in the slice (um)
slice_vector = 1:200;

% define the vector of times
time_vector = 1:100;

% allocate a vector to save the angles
angle_vector = zeros(size(slice_vector));

% allocate memory for the results
speed_map = zeros(length(time_vector),length(slice_vector));
% start counters for both variables
t_count = 1;
% for all the object distances
for time_v = time_vector
    a_count = 1;

    % for all the travel distances
    for travel = slice_vector
        % get the angle
        angle_vector(a_count) = travel/conversion_factor;
%         % get the distance travelled in visual space
%         visual_distance = 2*pi*distance*(angle_travel/360);
        % get the speed and save
        speed_map(t_count,a_count) = angle_vector(a_count)/(time_v*1e-3);
        % update the counter
        a_count = a_count + 1;
    end
    % update the counter
    t_count = t_count + 1;
end

% Plot
close all
figure

% define the density of labels
dens = 10;
% imagesc(log10(speed_map'))

% Rossi stimulus: 20deg/cyc, 2 cyc/s --> 40deg/s
rossi_speed = speed_map';
rossi_speed(round(rossi_speed)>39&round(rossi_speed)<41) = NaN;
imagesc(log10(rossi_speed))

set(gca,'TickLength',[0 0],'XTick',1:dens:length(time_vector),'XTickLabel',time_vector(1:dens:end),'XTickLabelRotation',0)
set(gca,'YTick',1:dens:length(slice_vector),'YTickLabel',slice_vector(1:dens:end))
xlabel('Displacement time (ms)')
ylabel('Slice distance (um)')

% superimpose the Rossi limits
hold on
min_dist = min_delta*conversion_factor;
plot(get(gca,'XLim'),[min_dist min_dist].*max(slice_vector)./length(slice_vector),'c')

max_dist = max_delta*conversion_factor;
plot(get(gca,'XLim'),[max_dist max_dist].*max(slice_vector)./length(slice_vector),'r')

plot([10 10].*max(time_vector)./length(time_vector), get(gca,'YLim'),'k')



set(gcf,'Color','w')
h = colorbar;
ylabel(h,'Speed log(deg/s)')
colormap(viridis)
%% Select the cells for comparison (from Simon)


od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
              ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci]]';
% a=[];
%read out DSI>0.25
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
%read out the actual direction for ipsi and contra
binodir=reshape([str(a).Dir],[2,length(a)])';
%get the delta 
binodir_delta=binodir(:,1)-binodir(:,2);
%find 'truly' binocular cells
idx_bi=find(abs(od_out_iviv(a,3))<0.1);
%remove cell if delta bino is more than 90 deg off
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90);
a(idx_bi(id_ov))=[];
%define direction sectors width
sector=60;
%define midpoint of widows
midpoint=130;
%define range
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
% g1=[];g2=[];g3=[];g4=[];
%find belonging to sectors 
g1=find(od_out_iviv(a,5)>s1a & od_out_iviv(a,5)<s1b) ;
g2=find(od_out_iviv(a,5)>s2a & od_out_iviv(a,5)<s2b);
g3=find(od_out_iviv(a,5)>s3a & od_out_iviv(a,5)<s3b);
g4=find(od_out_iviv(a,5)>s4a & od_out_iviv(a,5)<s4b);
%combine 315(M) and 135(L)
s1=[g1' g2'];
%combine 45(O) and 225(I)
s2=[g3' g4'];

% store the vectors
selection_both = {a(s1),a(s2)};

% get a flip vector (pref direction right and out)
flip_both = {zeros(length(s1),1),zeros(length(s2),1)};
flip_both{1}(1:length(g1)) = 1;
flip_both{2}(length(g2)+1:end) = 1;
%% Load a selection of maps
close all
% % define the target angles
% angle_vector = [45 225];
% angle_vector = [135 315];
% allocate memory for the traces
diff_trace = cell(2,1);
% for both directions
for dirs = 1:2
%     % define the target angle
%     target_angle = angle_vector(dirs);
%     angle_tolerance = 45;
%     
%     % generate a circular vector for DIRpref
%     dir_rotated = [str.DIRpref,[str.DIRpref]+360];
%     % get a selection vector
%     selection_vector = (dir_rotated > target_angle-angle_tolerance) & ...
%         (dir_rotated < target_angle+angle_tolerance);
%     % unrotate the selection vector
%     selection_vector = (selection_vector(1:70)+selection_vector(71:end))==1;
    selection_idx = selection_both{dirs};
    selection_vector = zeros(length(str),1);
    selection_vector(selection_idx) = 1;
    selection_vector = selection_vector==1;
    % get the maps
    maps_exc = cat(3,str(selection_vector).subpixel_excMap);
    maps_inh = cat(3,str(selection_vector).subpixel_inhMap);

    % get the number of maps
    number_maps = sum(selection_vector);
    %% Plot profiles

    figure

    % select the target layer
    target_layer = 2:6;
    % allocate memory for the deltas
    delta_vector = zeros(number_maps,2);
    % for all the maps
    for maps = 1:number_maps
        subplot(round(sqrt(number_maps)),ceil(sqrt(number_maps)),maps)
        
        % load the maps
        curr_exc = maps_exc(target_layer,:,maps);
        curr_inh = maps_inh(target_layer,:,maps);
        % if it's a flip, flip it
        if flip_both{dirs}(maps)
            curr_exc = fliplr(curr_exc);
            curr_inh = fliplr(curr_inh);
        end

%         exc_trace = [(sum(maps_exc(target_layer,:,maps),1)/sum(sum(maps_exc(target_layer,:,maps))))];
%         inh_trace = [(sum(maps_inh(target_layer,:,maps),1)/sum(sum(maps_inh(target_layer,:,maps))))];

        exc_trace = [0,diff(sum(curr_exc,1))];
        inh_trace = [0,diff(sum(curr_inh,1))];

        plot(exc_trace(1:8),'r')
        hold on
        plot(inh_trace(1:8),'b')
    %     plot(exc_trace(1:8)-inh_trace(1:8),'k')

        plot(-exc_trace(16:-1:9),'r--')
        plot(-inh_trace(16:-1:9),'b--')
%         plot((exc_trace(16:-1:9)-inh_trace(16:-1:9))-(exc_trace(1:8)-inh_trace(1:8)),'k--')
%         set(gca,'XTick',1:8,'XTickLabel',roundn(x_axis(1:8),-2),'XTickLabelRotation',45,'YLim',[-2 2])
        set(gca,'YLim',[-2 2])
        set(gcf,'Color','w')

        % save the delta
%         delta_vector(maps,:) = ...
%             (exc_trace(16:-1:9)-inh_trace(16:-1:9))-(exc_trace(1:8)-inh_trace(1:8));
%         delta_1 = exc_trace(1:8)-inh_trace(1:8);
%         delta_2 = -(exc_trace(16:-1:9)-inh_trace(16:-1:9));
%         [~,max_1] = max(abs(delta_1));
%         [~,max_2] = max(abs(delta_2));
%         delta_vector(maps,1) = max_1;
%         delta_vector(maps,2) = max_2;

        [~,max_exc_pref] = max(exc_trace(1:8));
        [~,max_inh_pref] = max(inh_trace(1:8));
        [~,max_exc_non] = min(exc_trace(16:-1:9));
        [~,max_inh_non] = min(inh_trace(16:-1:9));
        
        delta_vector(maps,1) = max_exc_pref - max_inh_pref;
        delta_vector(maps,2) = max_exc_non - max_inh_non;


        

    end
    % save the results for this angle
    diff_trace{dirs} = delta_vector;
    
end
%% Plot the traces
close all
figure

subplot(1,2,1)
plot(diff_trace{1}'.*69,'co-')
set(gca,'XLim',[0.5 2.5],'TickLength',[0 0],'YLim',[-2.1 2.1].*69,...
    'XTick',[1 2],'XTickLabel',{'Pref','Non-Pref'})
ylabel('Cortical distance (um)')

p1 = signrank(diff_trace{1}(:,1),diff_trace{1}(:,2),'tail','right');
title(strcat('Along slice p:',num2str(p1)))

% histogram(diff_trace{1}(:,1),'Normalization','Probability')
% hold on
% histogram(diff_trace{1}(:,2),'Normalization','Probability')

subplot(1,2,2)
plot(diff_trace{2}'.*69,'mo-')
set(gca,'XLim',[0.5 2.5],'TickLength',[0 0],'YLim',[-2.1 2.1].*69,...
    'XTick',[1 2],'XTickLabel',{'Pref','Non-Pref'})

p2 = signrank(diff_trace{2}(:,1),diff_trace{2}(:,2),'tail','left');
title(strcat('Across slice p:',num2str(p2)))

set(gcf,'Color','w')
%% Shuffle
% 
% % define the number of shuffles
% shuffle_num = 100;
% % allocate memory to hold the results
% shuffle_result = zeros(shuffle_num,2);
% % for all the shuffles
% for shuffle = 1:shuffle_num
% 
%     % allocate memory for the traces
%     diff_trace = cell(2,1);
%     % for both directions
%     for dirs = 1:2
%     %     % define the target angle
%     %     target_angle = angle_vector(dirs);
%     %     angle_tolerance = 45;
%     %     
%     %     % generate a circular vector for DIRpref
%     %     dir_rotated = [str.DIRpref,[str.DIRpref]+360];
%     %     % get a selection vector
%     %     selection_vector = (dir_rotated > target_angle-angle_tolerance) & ...
%     %         (dir_rotated < target_angle+angle_tolerance);
%     %     % unrotate the selection vector
%     %     selection_vector = (selection_vector(1:70)+selection_vector(71:end))==1;
% %         selection_idx = selection_both{dirs};
% 
%         % randomly pick a set of cells
%         selection_idx = randperm(length(str),30);
%         selection_vector = zeros(length(str),1);
%         selection_vector(selection_idx) = 1;
%         selection_vector = selection_vector==1;
%         % get the maps
%         maps_exc = cat(3,str(selection_vector).subpixel_raw_excMap);
%         maps_inh = cat(3,str(selection_vector).subpixel_raw_inhMap);
% 
%         % get the number of maps
%         number_maps = sum(selection_vector);
% 
%         % select the target layer
%         target_layer = 2:6;
%         % allocate memory for the deltas
%         delta_vector = zeros(number_maps,2);
%         % for all the maps
%         for maps = 1:number_maps
% 
%             exc_trace = [0,diff(sum(maps_exc(target_layer,:,maps),1)/sum(sum(maps_exc(target_layer,:,maps))))];
%             inh_trace = [0,diff(sum(maps_inh(target_layer,:,maps),1)/sum(sum(maps_inh(target_layer,:,maps))))];
% 
% 
% 
%             % save the delta
%             delta_vector(maps,1) = max(abs(exc_trace(1:8)-inh_trace(1:8)));
%             delta_vector(maps,2) = max(abs(exc_trace(16:-1:9)-inh_trace(16:-1:9)));
%         end
%         % save the results for this angle
%         diff_trace{dirs} = delta_vector;
% 
%     end
%     
%     % store the shuffle result
%     shuffle_result(shuffle,:) = [signrank(diff_trace{1}(:,1),diff_trace{1}(:,2),'tail','left'),...
%         signrank(diff_trace{2}(:,1),diff_trace{2}(:,2),'tail','left')];
% end
% %% Plot the shuffles
% close all
% figure
% 
% subplot(1,2,1)
% % plot(diff_trace{1}','co-')
% histogram(-log10(shuffle_result(:,1)),20,'Normalization','Probability');
% hold on
% plot([1.3 1.3],get(gca,'YLim'),'r')
% set(gca,'TickLength',[0 0])
% panel1 = gca;
% ylabel('Counts')
% xlabel('-log10(p value)')
% title('Along slice')
% 
% % histogram(diff_trace{1}(:,1),'Normalization','Probability')
% % hold on
% % histogram(diff_trace{1}(:,2),'Normalization','Probability')
% 
% subplot(1,2,2)
% 
% histogram(-log10(shuffle_result(:,2)),20,'Normalization','Probability');
% hold on
% set(gca,'TickLength',[0 0],'YLim',get(panel1,'YLim'))
% 
% plot([1.3 1.3],get(gca,'YLim'),'r')
% xlabel('-log10(p value)')
% title('Across slice')
% 
% set(gcf,'Color','w')
% % histogram(diff_trace{2}(:,1),'Normalization','Probability')
% % hold on
% % histogram(diff_trace{2}(:,2),'Normalization','Probability')
% 
% % subplot(1,2,1)
% % plot(diff_trace{1}','c--')
% % hold on
% % plot(nanmean(diff_trace{1},1),'c', 'LineWidth',2)
% % set(gca,'XTick',1:8,'XTickLabel',roundn(x_axis,-2),'XTickLabelRotation',45)
% % set(gca,'YLim',[-0.2 0.16])
% % subplot(1,2,2)
% % plot(diff_trace{2}','m--')
% % hold on
% % plot(nanmean(diff_trace{2},1),'m', 'LineWidth',2)
% % set(gca,'YLim',[-0.2 0.16])
% % set(gca,'XTick',1:8,'XTickLabel',roundn(x_axis,-2),'XTickLabelRotation',45)



