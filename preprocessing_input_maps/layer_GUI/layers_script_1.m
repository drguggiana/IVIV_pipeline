clearvars
close all
% get the paths
Paths
%select the folders to process
folder_list = uipickfiles('FilterSpec',input_maps_path);

%define the path to save the output files
out_path = layers_path;
%define the path to the default grid
default_grid = default_grid_path;
%get the number of folders
folder_num = length(folder_list);
%allocate memory to store the subfolders that qualify
folder_cell = cell(folder_num,1);
%allocate a vector to hold the pixel size of the image based on the setup
%(x,y)
pixel_size = zeros(folder_num,2);
%also allocate to store the pixels sizes for each subfolder
pixel_cell = cell(folder_num,1);

%for all the folders
for folders = 1:folder_num

    %identify the setup for this folder
    setup_id = strsplit(folder_list{folders},'\');
    setup_id = str2double(setup_id{end-1}(end));
    %define the image format and pixel size in microns based on the setup
    switch setup_id
        case 1
            im_format = '*.bmp';
            pixel_size(folders,:) = [1580,1183];
        case 2
            im_format = '*.tif';
            pixel_size(folders,:) = [2548,1903];
    end
    
    %get the subfolders in the first level of this folder
    subfolder_list = dir(folder_list{folders});

    %also, get rid of the non-folders
    subfolder_list = subfolder_list(vertcat(subfolder_list(:).isdir)==1);
    %leave only the folders that have SW or FM in them
    subfolder_list = subfolder_list(contains({subfolder_list(:).name},'SW')|contains({subfolder_list(:).name},'FM'));
    %finally, build the full paths for the folders, including the "images"
    %folder and the "map01"
    path_list = cell(size(subfolder_list,1),2);
    %create a vector to signal which fields to eliminate later
    elim_vec = ones(size(path_list,1),1);
    %for all the subfolders
    for subs = 1:size(subfolder_list,1)
        path_list{subs,1} = strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\images\');
        temp_list = dir(strcat(path_list{subs,1},im_format));
        if isempty(temp_list)
            elim_vec(subs) = 0;
            continue
        end
        path_list{subs,1} = strcat(path_list{subs,1},temp_list(1).name);
        path_list{subs,2} = strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\map*\');
        temp_list = dir(strcat(path_list{subs,2},'*.xsg'));
        path_list{subs,2} = strcat(temp_list(1).folder,'\',temp_list(1).name);
    end
    %update the folder cell only with the cells that had an image
    folder_cell{folders} = path_list(elim_vec==1,:);
    %and fill in the pixel sizes accordingly
    pixel_cell{folders} = ones(length(path_list(elim_vec==1,:)),2).*pixel_size(folders,:);
end

%concatenate the entire list of paths
folder_all = vertcat(folder_cell{:});
%and the pixel sizes
pixel_all = vertcat(pixel_cell{:});
%get the number of cells to do
cell_num = size(folder_all,1);
%run the GUI for all the cells, save in every iteration
%for all the cells
for cells = 1:cell_num
    %run the GUI with the paths from this cell
    gui_handle = layer_GUI_1([],folder_all{cells,1},folder_all{cells,2},cells,cell_num,pixel_all(cells,:),default_grid);
    waitfor(gui_handle)
    %save the output in a unique file for this cell
    %get the cell's name
    cell_name = strsplit(folder_all{cells,1},'\');
    cell_name = strcat(cell_name{6}(1:2),cell_name{5}(1:6),cell_name{6}(3:6));
    %save the variable
    save(fullfile(out_path,cell_name),'new_grid')
    
end