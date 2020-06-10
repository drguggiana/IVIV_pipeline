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


%% Load python (only needed if it changes)
% pyversion(python_exe_path)
% pyversion
 
insert(py.sys.path,int32(0),'C:\Users\drguggiana\PycharmProjects\L23_Pipeline\code');
%% Get the data to be clustered and the linkage matrix
close all

% get the fractions
frv = reshape([str(:).frac_vert],32,length(str))';
L23fr = [sum(frv(:,3:5),2) sum(frv(:,19:21),2)];
L4fr = [sum(frv(:,6:7),2) sum(frv(:,22:23),2)];
L5fr = [sum(frv(:,8:11),2) sum(frv(:,24:27),2)];

input_data = [L23fr L4fr L5fr];


% get the pia_input
pia_input = [str(:).pialD]';
% define the number of clusters

clu_num =2;
% cluster
[idx_input, clustering_input, leafOrder] = hca([L23fr L4fr L5fr]...
    ,0,'ward',clu_num,pia_input,0,0.6);
%% Call the node significance from python

py.dendrogram_significance

clu_num = 8;

% linkage method list
link_list = {'ward','complete','average'};
% allocate memory for the indexes
index_cell = cell(length(link_list),1);

% for all linkage methods
for linkage_method = 1
    % cluster
    [index_cell{linkage_method}, clustering_input, leafOrder] = hca(input_data...
        ,0,link_list{linkage_method},clu_num,pia_input,0,0.6);

    %% Import the module

    mod = py.importlib.import_module('dendrogram_significance');
    %% Run the function

    results = cell(mod.dendrogram_significance(input_data(:)',...
        pyargs('alpha',0.05,'link',link_list{linkage_method})));


    distance_result = cellfun(@(x) str2double(char(x)),cell(results{1}),'UniformOutput',false);
    distance_result = vertcat(distance_result{:});

    ci_result = cellfun(@(x) str2double(char(x)),cell(results{1}),'UniformOutput',false);
    ci_result = vertcat(ci_result{:});
    %% Plot the dendrogram
%     close all
    figure
    [H,T,outperm] = dendrogram(clustering_input,200);
    % recolor the nodes
    % for all the elements
    for el = 1:size(H,1)
        if distance_result(el) == 0
            H(el).Color = [0 0 0];
        else
            H(el).Color = [1 0 0];
        end
    end
end
%% Get the indexes for kmeans

[idx_k,C,SUMD,K,PC] = kmeans_opt([L23fr L4fr L5fr],10,0.95,1000); 
%% Quantify the overlap between index lists

% define the idx to use
idx1 = index_cell{1};
idx2 = idx_k;
% get the cluster numbers
clu_num1 = length(unique(idx1));
clu_num2 = length(unique(idx2) );

% define the noise threshold (in percentage from the total number)
noise_threshold = 1/length(idx1);
% calculate the overlap
[overlap,pairs,ov_matrix] = quantify_overlap(idx1,idx2,noise_threshold);
%% compare to shuffle
% close all
% define the number of shuffles
shuffle_number = 1000;

% allocate memory for the results
shuffle_results = zeros(shuffle_number,1);
% for all the shuffles
for shuff = 1:shuffle_number
    % shuffle the idx2
%     random_idx1 = idx1(randperm(length(idx1)));
    random_idx2 = idx2(randperm(length(idx2)));
    % calculate the overlap
    [shuffle_results(shuff),~,a] = quantify_overlap(idx1,random_idx2,noise_threshold);
end

% plot the distribution
figure
histogram(shuffle_results,'Normalization','Probability')
hold on
plot([prctile(shuffle_results,5),prctile(shuffle_results,5)],get(gca,'YLim'),'-k')
plot([prctile(shuffle_results,95),prctile(shuffle_results,95)],get(gca,'YLim'),'-k')
plot([overlap,overlap],get(gca,'YLim'),'-r')
%% Plot the optimal pairs
% close all
figure
% sort the rows
idx_sort = sortrows([idx1,idx2]);
% plot
subplot(2,1,1)
imagesc(idx_sort(:,1)')
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])
ylabel('HCA')
cba = colorbar;
set(cba,'TickLength',0,'YTick',1:length(unique(idx1)))
colormap(magma)

subplot(2,1,2)
imagesc(idx_sort(:,2)')
set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])
ylabel('K-means')
cba = colorbar;
set(cba,'TickLength',0,'YTick',1:length(unique(idx2)))
colormap(magma)

set(gcf,'Color','w')
% suptitle('Cluster Overlap')
% imagesc(1:length(idx1))

