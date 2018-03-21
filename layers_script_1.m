clearvars
close all

%select the folders to process
folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\EXPLORER ONE');

%define the path to save the output files
out_path = 'R:\Share\Simon\Drago_Volker_Simon\layer_GUI_out';
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
    path_list = cell(size(subfolder_list,1),2);
    %create a vector to signal which fields to eliminate later
    elim_vec = ones(size(path_list,1),1);
    %for all the subfolders
    for subs = 1:size(subfolder_list,1)
        path_list{subs,1} = strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\images\');
        temp_list = dir(strcat(path_list{subs,1},'*.tif'));
        if isempty(temp_list)
            elim_vec(subs) = 0;
            continue
        end
        path_list{subs,1} = strcat(path_list{subs,1},temp_list(1).name);
        path_list{subs,2} = strcat(subfolder_list(subs).folder,'\',subfolder_list(subs).name,'\map01\');
        temp_list = dir(strcat(path_list{subs,2},'*.xsg'));
        path_list{subs,2} = strcat(path_list{subs,2},temp_list(1).name);
    end
    %update the folder cell only with the cells that had an image
    folder_cell{folders} = path_list(elim_vec==1,:);
end

%concatenate the entire list of paths
folder_all = vertcat(folder_cell{:});
%get the number of cells to do
cell_num = size(folder_all,1);
%run the GUI for all the cells, save in every iteration
%for all the cells
for cells = 1:cell_num
    %run the GUI with the paths from this cell
    gui_handle = layer_GUI_1([],folder_all{cells,1},folder_all{cells,2},cells,cell_num);
    waitfor(gui_handle)
    %save the output in a unique file for this cell
    %get the cell's name
    cell_name = strsplit(folder_all{cells,1},'\');
    cell_name = strcat(cell_name{4},'_',cell_name{5});
    %save the variable
    save(fullfile(out_path,cell_name),'new_grid')
    
end