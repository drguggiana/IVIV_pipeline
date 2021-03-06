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

% get the pial depth
pialD=cat(1,str.pialD);

% get the layer 4 excitation
frac4ex = cat(1,str.frac_vert);
frac4ex = sum(frac4ex(:,6:7),2);

% get the inhibitory L23 angle
out_ang_inL23 = cat(1,str.ang_inL23);
% get the inhibitory L23 x centroid
centroidX23in=abs(out_ang_inL23(:,3)-out_ang_inL23(:,1));
centroidY23in=abs(out_ang_inL23(:,4)-out_ang_inL23(:,2));

% assemble the feature vector
cell_cell = cat(2,pialD,frac4ex,centroidY23in,centroidX23in);

% cell_cell = cat(2,pcs(:,2),ang);
cell_cell = normr_2(cell_cell,2);
%% Run UMAP on the data (or load the embedding)

% run the embedding from scratch
% load the embedding
reduced_data = load(umap_path);
reduced_data = reduced_data.reduced_data;
% [reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 15, 'min_dist', 0.5);

% % load the umap file from path
% reduced_data = load(umap_path);
% reduced_data = reduced_data.reduced_data;
%% Plot the umaps

close all

% include the field names that should be plotted. After the plot list, the
% other inputs determine which fields are log scaled and normalized
% respectively (vector)
% plot_list = {'hemisphere','sliceOri','pialD','Cluster_id','somaCenter','cellID'};

% plot_list = {'morph'};

% plot_list = {'OSIpref','DSIpref','ODIpref','Capeakpref','ORIpref','DIRpref','noise','PCs',...
%     'ang_exL4','ang_inL4','pialD','pci',...
%     'frac_vert','corr_exc_apical','corr_exc_basal','corr_inh_apical','corr_inh_basal'};
plot_list = {'OSIpref','DSIpref','ODIpref','Sigmapref','Capeakpref','ORIpref','DIRpref','noise','PCs',...
   'ang_exL23','ang_inL23','pialD','pci',...
   'frac_vert','resp'};

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
