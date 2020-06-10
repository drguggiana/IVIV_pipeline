
%% 
clu_num =2;
%pcs =[];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input_ward, clustering_input, leafOrder] = hca([L23fr L4fr L5fr],0,'ward',clu_num,pia_input,3,0.6);%call function for clustering
%[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%[idx_input_comp, clustering_input, leafOrder] = hca([L23fr L4fr L5fr],0,'complete',clu_num,pia_input,3,0.6);%call function for clustering
%% Plot average maps of clusters
plot_avg_maps(str,1:147,ex_map,in_map,pia_input,10,1,idx_input);
%% 
%% Looking at input fractions
dendroplot_SW(clustering_input,leafOrder,2,[L23fr(:,1) L4fr(:,1)  L5fr(:,1)],{'L2/3ex','L4ex','L5ex','L2/3in','L4in','L5in'})
%% 
dendroplot_SW(clustering_input,leafOrder,2,[L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3ex','L4ex','L5ex','L2/3in','L4in','L5in'},pia_input)
%% kmeans clustering
idx_input=[];
 idx_input2 = kmeans([L23fr L4fr L5fr],5)
 %% 
 
 %% 
 
 [idx_k,C,SUMD,K,PC]=kmeans_opt([L23fr L4fr L5fr],50);

%% Barplot difference of clusters
%Pial depth
[statsout] = barplot_sw(pia_input,idx_input,{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45)
%% 
[statsout] = barplot_sw(span(:,1),idx_input,{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45)
%% 
[statsout] = barplot_sw(abs(ang_inL23(:,4)-ang_inL23(:,2)),idx_input,{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45)
%% 

a=find(od_out_iviv(:,1)>0.25)  
[statsout] = barplot_sw(od_out_iviv(a,4),idx_input(a),{'Clusters','Orientation preference'});xtickangle(45);set(gca,'Ydir','reverse');set(gca,'FontSize',12)
%% 
a=find(od_out_iviv(:,2)>0.25);  
par=od_out_iviv(a,9)
g1=find(idx_input(a)==1) ;
g2=find(idx_input(a)==2);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Pial depth (µm)','Color','k'); ;set(gca,'FontSize',10);