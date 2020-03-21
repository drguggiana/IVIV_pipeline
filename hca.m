function [idx, clustering, leafOrder] = hca(input,pc_ana,linkage_meth,cluster_number,pia,inputd,threshold)
if inputd==1;
Z=normr_2(input,2);
else
 Z=zscore(input);
end
if pc_ana==1
[coeff,score,latent,tsquared,explained,mu] = pca(Z);

nr_PCs=find(cumsum(explained)>80,1);
score_s=score(:,1:nr_PCs);
data_ve=score_s;
D=pdist(score_s);%Euclidean distance between pairs of observations
clustering=linkage(D,linkage_meth);%complete:Farthest distance %ward:Inner squared distance (minimum variance algorithm)
leafOrder = optimalleaforder(clustering,D);%returns an optimal leaf ordering for the hierarchical binary cluster tree
idx=cluster(clustering,'maxclust',cluster_number);
else
data_ve=Z;
D=pdist(Z);%Euclidean distance between pairs of observations
clustering=linkage(D,linkage_meth);%complete:Farthest distance %ward:Inner squared distance (minimum variance algorithm)
leafOrder = optimalleaforder(clustering,D);%returns an optimal leaf ordering for the hierarchical binary cluster tree
idx=cluster(clustering,'maxclust',cluster_number);
end

% figure();
% set(gcf,'color','w');
%  %str={'V_{min}(mV)','V_{peak}(mV)','V_{thresh}(mV)', 'Vslope_{max} (\DeltamV/\Deltams)','V_{half} (mV)','Spike_{amplitude} (mV)','AHP_{max}(mV)','Spike_{rise} (ms)', 'Spike_{fall} (ms)','Spike_{base width}(ms)', 'Spike_{half width} (ms)','First spike latency (ms)', 'V_{rest} (mV)','Tau (ms)','R_{IN}(M\Omega)','Sag ratio (%)','Rheobase (pA)','Spike frequency_{max} (Hz)','Pial depth (µm)'}
% ax=max(input);
% axm=min(input);
% for i=1:size(input,2);
% hold on;
% subplot(4,5,i)
% h=histogram(input(:,i),'LineWidth',1.5,'DisplayStyle','stairs')
% h.EdgeColor = 'k';
% %h.FaceColor = [0.5 0.5 0.5];
% set(gca,'box','off');
% %xlabel(str(i));
% %ylim([0 80]);
% if ax(i)<0
% xlim([axm(i)+ axm(i)*0.25 0]);
% else
%     xlim([0 ax(i)+ ax(i)*0.25]);
% 
% end
% end

color_idx={'m','g','b','r'};
fig2 = figure;
set(fig2, 'Name', 'HCA');
set(fig2, 'Position', [200, 0, 1500, 1000]);
set(gcf,'color','w');

subplot(3,2,1);
H=dendrogram(clustering,0,'reorder',leafOrder,'Orientation','left','ColorThreshold',[threshold*max(clustering(:,3))]);%plot entire dendrogram (0), with the reordered leafOrder
%dendrogram(clustering,0,'reorder',leafOrder)
set(gcf,'color','w');
set(H,'LineWidth',1.5);
% ylabel('Euclidean distance');
% xlabel('Cells');
 xlabel('Euclidean distance');
 %ylabel('Cells');
set(gca,'box','off');
%axis square;
 lineColours = cell2mat(get(H,'Color'));
%set(gca, 'XTick',[]);
set(gca, 'YTick',[]);
% H(1).Color = 'r';
% H(8).Color = [0 0.5 0];
for m=1:length(cluster_number)
idx_c=find(lineColours(:,m)==1);
for i=1:length(idx_c)
    H(idx_c(i)).Color=color_idx{m};
end
idx_c=[];
% idx2=find(lineColours(:,2)==1);
% idx1=find(lineColours(:,3)==1);
end

% 
% 
% for i=1:length(idx2)
%     H(idx2(i)).Color=[1 1 0.1];
% end
% 
% for i=1:length(idx3)
%     H(idx3(i)).Color=[0 1 0];
% end




%Siluohtte test for HCA, cluster number confirmation
subplot(3,2,2);
for i=1:10;  
idx_HCA = cluster(clustering,'maxclust',i);
sil=silhouette(data_ve,idx_HCA);
avg_sil=mean(sil);
plot(i,avg_sil,'m--o');
set(gca,'box','off');
axis square;
xlabel('Number of clusters');
ylabel('Average silhouette value');
hold on;
end

%Thorndike procedure for HCA, cluster number confirmation
subplot(3,2,3);plot([2:10],flip(diff(clustering(end-9:end,3))),'-o','Color','m'); ylabel('\Delta Euclidean distance between clusters joined');
set(gcf,'color','w');
xlabel('Number of clusters');
set(gca,'box','off');
axis square;xlim([0 10]);

%Elbow procedure for Kmeans, cluster number confirmation
[idx_k,C,SUMD,K,PC]=kmeans_opt(data_ve,10);
subplot(3,2,4);
plot(PC,'-o','Color','m');
set(gcf,'color','w');
xlabel('Number of clusters');
ylabel('% variance');
set(gca,'box','off');
axis square;

idxn=[min(idx):max(idx)];


subplot(3,2,5);
hold on;
for m=1:length(idxn)
clus_pia{:,m}=pia(find(pia(idx==m)));
%barwitherr(m,nanstd(clus_pia{:,m})/sqrt(length(clus_pia{:,m})),nanmean(clus_pia{:,m}));
 errorbar(m,nanmean(clus_pia{:,m}),nanstd(clus_pia{:,m})/sqrt(length(clus_pia{:,m})));
hold on;
xlim([0 max(length(idxn))+1]);
ylim([min(clus_pia{:,m}) max(clus_pia{:,m})]);
axis square;
set(gca,'box','off');
set(gcf,'color','w');
ylabel('pial depth');
xlabel('cluster');
%hold on;scatter(repmat(m,length(clus_pia{:,m}),1),clus_pia{:,m});
end



end