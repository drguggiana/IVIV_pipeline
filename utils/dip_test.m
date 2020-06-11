%% Hartigan's dip test with structure data
%% load the paths and clean up
clearvars
close all

Paths
%% Do a Hartigan's Dip test
close all

% define the number of bootstraps
nboot = 500;
% define the data to be tested
 [coeff,score,latent,~,explained,mu] = pca([L23fr(:,1) L4fr(:,1) L5fr(:,1) L23fr(:,2) L4fr(:,2) L5fr(:,2)]);
 var_exp(explained,[],[]); 
 ylim([0 100]);;xlim([0 7]);xticks([1:1:6]);
 yticks([0:25:100])
%  
%% 
test_data =score(:,1:5);
%test_data = vertcat(str.PCs);
% get the number of dimensions
data_dim = size(test_data,2);

% allocate memory for the results
dip_results = zeros(data_dim,2);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 600, 300]);
% check for all PCs
for i = 1:data_dim
    [dip_results(i,1), dip_results(i,2)] = HartigansDipSignifTest(test_data(:,i), nboot);
    
    % plot the result
    subplot(round(sqrt(data_dim)),ceil(sqrt(data_dim)),i)
   h1= histogram(test_data(:,i));box off;h1.FaceColor='w';
   ylim([0 70]);yticks([0:35:70])
    title(strjoin({'Dip test:',num2str(dip_results(i,1)),'p val:',num2str(dip_results(i,2))},' '))
end