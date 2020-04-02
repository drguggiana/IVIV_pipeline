% produce UMAP plots of the data

%% load the paths and clean up
clearvars
close all

Paths
%% Load the relevant files

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Get rid of the entries without PCs

% get the field contents as a cell
pc_vector = {str.PCs};
% turn into a logical vector
pc_vector = ~cellfun(@isempty,pc_vector);
% kill the fields
str = str(pc_vector);
%% Get the PCs and assemble the feature vector

pcs = cat(1,str.PCs);
cell_cell = cat(2,pcs(:,1:3),pcs(:,4:6));
%% Run UMAP on the data

[reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 10, 'min_dist', 0.1);
%% Plot the umaps

close all

% include the field names that should be plotted. After the plot list, the
% other inputs determine which fields are log scaled and normalized
% respectively (vector)
plot_list = {'OSIpref','DSIpref','ODIpref','ORIpref','DIRpref','noise','PCs',...
    'Cluster_id','ang_exL23','Capeakpref','pci','sad'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula')

autoArrangeFigures


