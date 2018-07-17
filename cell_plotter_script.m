%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load the data structure

%define the main file path
model_path = 'R:\Share\Simon\Drago_Volker_Simon\Full_data_structure_joel';

%load the file
%define the file tags
file_tags = {'_dataStruct'};
%allocate memory to store the data
all_data = cell(length(file_tags),1);

%for all the tags
for tags = 1:length(file_tags)
    %load the most recent file from the overlap ones
    listing = dir(strcat(model_path,'\*',file_tags{tags},'*.mat'));
    dates = datetime({listing.date});
    [~,ind] = max(dates);
    load_name = listing(ind).name;
    
    %load all the variables in the file
    load(fullfile(model_path,load_name));
end
%get the number of cells
cell_num = length(invitro_struct);
%% Plot a subset of cells

close all
%define the number of cells to plot
cell_toplot = 1;

%get the indexes
%cell_idx = randperm(cell_num,cell_toplot);
cell_idx = 1;
%define the smoothing factor
sf = 1;

%for all the cells to plot
for cells = 1:cell_toplot
    
    %plot the excitatory map
    %get the map
    exc_map = invitro_struct(cell_idx(cells)).excMap;
    
    %check if the map exists, otherwise skip
    if ~isnan(sum(exc_map(:)))
        %create the figure to put the plot on
        h = figure;
        %define the plot type (2 for excitatory)
        plot_type = 2;
        %get the pial distance
        soma_info = invitro_struct(cell_idx(cells)).somaCenter;
        %plot the map
        map_plot3(exc_map,'',plot_type,h,sf,0,1,soma_info);
    end
    
    %plot the inhibitory map
    %get the map
    inh_map = invitro_struct(cell_idx(cells)).inhMap;
    
    %check if the map exists, otherwise skip
    if ~isnan(sum(inh_map(:)))
        %create the figure to put the plot on
        h2 = figure;
        %define the plot type (3 for inhibitory)
        plot_type = 3;
        %get the pial distance
        soma_info = invitro_struct(cell_idx(cells)).somaCenter;
        %plot the map
        map_plot3(inh_map,'',plot_type,h2,sf,0,1,soma_info);
    end
    
    %plot the overlap map
    %get the map
    ove_map = cat(3,invitro_struct(cell_idx(cells)).excMap,invitro_struct(cell_idx(cells)).inhMap);
    
    %check if both maps exist, otherwise skip
    if ~isnan(sum(exc_map(:)))&&~isnan(sum(inh_map(:)))
        %create the figure to put the plot on
        h3 = figure;
        %define the plot type (1 for overlap)
        plot_type = 1;
        %get the pial distance
        soma_info = invitro_struct(cell_idx(cells)).somaCenter;
        %plot the map
        map_plot3(ove_map,'',plot_type,h3,sf,0,1,soma_info);
    end
end
%% OFF Compare soma info and pialD variables
% close all
% 
% for cells = 1:cell_num
%     plot(invitro_struct(cells).somaCenter(2),invitro_struct(cells).pialD,'*')
%     hold('on')
%     
% end
% 
% xlabel('Soma Center')
% ylabel('pialD')
% 
% somacent = cat(1,invitro_struct(:).somaCenter);
% somacent = somacent(:,2);
% pial = cat(1,invitro_struct(:).pialD);
% 
% corr(somacent,pial)
% 
% mean(somacent-pial)
% std(somacent-pial)
% figure
% histogram(somacent-pial)