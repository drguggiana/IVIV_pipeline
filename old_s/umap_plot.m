function umap_plot(data, vis_data, desc)

%% UMAP embedding



cell_cell = data;
% cell_cell = cat(2,pcs(:,2),ang);
cell_cell = normr_2(cell_cell,2);
%% Run UMAP on the data
% run the embedding from scratch
[reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 30, 'min_dist', 0.5);


%% 
if size(vis_data,1)==1
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 350, 300]);
set(gcf,'color','w');
hold on
scatter(cell_cell(:,1),cell_cell(:,2),30,vis_data,'filled')
;hold on;colorbar;xlabel('dim1');ylabel('dim2');
title(desc);
else
    
end
end