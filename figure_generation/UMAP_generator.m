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
% cell_cell = cat(2,pcs(:,1:3),pcs(:,4:6));
cell_cell = pcs;
frac = cat(1,str.frac_vert);
frac4ex = mean(frac(:,6:7),2);
frac4in = mean(frac(:,22:23),2);
frac23ex = mean(frac(:,3:5),2);
frac23in = mean(frac(:,19:21),2);

ang23ex = cat(1,str.ang_exL23);
alph23ex = ang23ex(:,5);
bet23ex = ang23ex(:,6);
centroidX23ex = ang23ex(:,3);

ang23in = cat(1,str.ang_inL23);
alph23in = ang23in(:,5);
centroidX23in = ang23in(:,3);

soma = cat(1,str.subpixel_soma);
pialD = cat(1,str.pialD);
% cell_cell = cat(2,pialD,frac4ex,bet23ex,alph23ex);
cell_cell = cat(2,pialD,frac23ex,frac23in,frac4ex,frac4in,alph23ex,alph23in,...
    centroidX23ex,centroidX23in);

% cell_cell = cat(2,pcs(:,2),ang);
cell_cell = normr_2(cell_cell,2);
%% Run UMAP on the data

[reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 15, 'min_dist', 0.5);
%% Plot the umaps

close all

% include the field names that should be plotted. After the plot list, the
% other inputs determine which fields are log scaled and normalized
% respectively (vector)
% plot_list = {'hemisphere','sliceOri','pialD','Cluster_id','somaCenter','cellID'};

plot_list = {'OSIpref','DSIpref','ODIpref','ORIpref','DIRpref','noise','PCs',...
    'Cluster_id','ang_exL23','ang_inL23','sliceOri','pialD','pci',...
    'frac_vert','corr_exc_apical','corr_exc_basal','corr_inh_apical','corr_inh_basal','subpixel_soma'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula')

autoArrangeFigures
%%
error('stop here')
%%
close all
soma = cat(1,str.subpixel_soma);
soma = soma(:,1);
slice_ori = cat(1,str.sliceOri);
cluster_idx = cat(1,str.Cluster_id);
ang = cat(1,str.ang_exL23);
ang = ang(:,3);

figure
% histogram(soma(slice_ori==1),10)
% hold on
% histogram(soma(slice_ori==0),10)

% plot(soma(cluster_idx==1),slice_ori(cluster_idx==1),'o')
% mean(soma(cluster_idx==1))
% hold on
% plot(soma(cluster_idx==3),slice_ori(cluster_idx==3),'o')
% mean(soma(cluster_idx==3))

cell_vec = (1:147)';
plot(soma,ang,'o')
% plot(soma(cluster_idx==1&cell_vec<46),ang(cluster_idx==1&cell_vec<46),'bo')
% hold on
% plot(soma(cluster_idx==1&cell_vec>=46),ang(cluster_idx==1&cell_vec>=46),'b*')
% plot(-soma(cluster_idx==2&cell_vec<46),ang(cluster_idx==2&cell_vec<46),'ro')
% plot(-soma(cluster_idx==2&cell_vec>=46),ang(cluster_idx==2&cell_vec>=46),'r*')
xlabel('soma x')
ylabel('centroid x')
axis square
%% 
% close all
var1 = soma(:,1);
var2 = pcs(:,2);
figure
plot(var1,var2,'o')

non_nan = ~(isnan(var1)|isnan(var2));
[rho,pval] = corr(var1(non_nan),var2(non_nan));
title(strjoin({num2str(rho),num2str(pval)},'_'),'Interpreter','None')
axis square
