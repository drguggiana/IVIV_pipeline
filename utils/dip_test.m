%% Hartigan's dip test with structure data
%% load the paths and clean up
clearvars
close all

Paths

% load the main structure
main_path = structure_file_path;
str = load(main_path);
str = str.str;
%% Do a Hartigan's Dip test
close all

% define the number of bootstraps
nboot = 500;
% define the data to be tested
% test_data = input_data;
test_data = vertcat(str.PCs);
% get the number of dimensions
data_dim = size(test_data,2);

% allocate memory for the results
dip_results = zeros(data_dim,2);
figure
% check for all PCs
for i = 1:data_dim
    [dip_results(i,1), dip_results(i,2)] = HartigansDipSignifTest(test_data(:,i), nboot);
    
    % plot the result
    subplot(round(sqrt(data_dim)),ceil(sqrt(data_dim)),i)
    histogram(test_data(:,i))
    title(strjoin({'Dip test:',num2str(dip_results(i,1)),'p val:',num2str(dip_results(i,2))},' '))
end