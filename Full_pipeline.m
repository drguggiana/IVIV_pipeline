%% Main script to run the entire pipeline (just press F5 and select the files, first Setup1 then Setup2)

%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load paths

%allocate memory to store the paths and the ori_setup from both setups
path_info = cell(2,2);

%for both setups
for setups = 1:2

    %pick the folders to use
    folder_list = uipickfiles('FilterSpec','I:\Simon Weiler\INPUT MAPS_final\');


    %define the path to save the output files
    out_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_cluster_out';
    %get the number of folders
    folder_num = length(folder_list);
    %allocate memory to store the subfolders that qualify
    folder_cell = cell(folder_num,1);
    %get the setup of origin (for later)
    ori_setup = strsplit(folder_list{1},'\');
    ori_setup = ori_setup{end-1};

    %for all the folders
    for folders = 1:folder_num

        %get the subfolders in the first level of this folder
        subfolder_list = dir(folder_list{folders});

        %also, get rid of the non-folders
        subfolder_list = subfolder_list(vertcat(subfolder_list(:).isdir)==1);
        %leave only the folders that have SW in them
        subfolder_list = subfolder_list(contains({subfolder_list(:).name},{'SW','FM'}));
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
    
    %store the in the info cell
    path_info{setups,1} = folder_all;
    path_info{setups,2} = ori_setup;
end
%% Run the interpolation software
%for both setups
for setups = 1:2
    PIPE_Trace_interpolate(path_info{setups,1},path_info{setups,2})
end
%% Run the mapmaker software

%define the vector of file tags
file_tags = {{'_rawClassSetup1'},{'_rawClassSetup2'}};
%for both setups
for setups = 1:2
    PIPE_Trace_mapmaker_Class(file_tags{setups})
end
%% Run the structure creator software

PIPE_Map_structureCreator