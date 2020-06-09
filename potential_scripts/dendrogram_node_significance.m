%% Perform significance node analysis
%% load the paths and clean up
clearvars
close all

Paths
%% Load the data
% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Get the data to be clustered and the linkage matrix

% get the fractions
frv = reshape([str(:).frac_vert],32,length(str))';
L23fr = [sum(frv(:,3:5),2) sum(frv(:,19:21),2)];
L4fr = [sum(frv(:,6:7),2) sum(frv(:,22:23),2)];
L5fr = [sum(frv(:,8:11),2) sum(frv(:,24:27),2)];

% get the pia_input
pia_input = [str(:).pialD]';
% define the number of clusters
clu_num =2;
% cluster
[idx_input, clustering_input, leafOrder] = hca([L23fr L4fr L5fr]...
    ,0,'ward',clu_num,pia_input,0,0.6);
%% Call the node significance from python

py.dendrogram_significance