%% Clean up ZAZZZZZZZ
clearvars
close all
%% Load files

%define the file tags
file_tags = {'_invitro','_invivo','_morpho'};

%for all the tags
for tags = 1:length(file_tags)
    %loading path (might have to edit this, depending on your own path)
    load_path = 'R:\Share\Simon\Drago_Volker_Simon\Analysis files';
    %load the most recent file from the overlap ones
    listing = dir(strcat(load_path,'\*',file_tags{tags},'.mat'));
    dates = datetime({listing.date});
    [~,ind] = max(dates);
    load_name = listing(ind).name;
    
    %load all the variables in the file
    load(fullfile(load_path,load_name))
end
%% OFF Compare the connections between the different clustering systems

% close all
% 
% %invitro and morpho
% 
% [~,ia,ib] = intersect(invitro_names,morpho_names);
% 
% %create averages of the common cell import properties using the local
% %morpho_clusters
% %get the number of common cells
% com_num = length(ia);
% %get a vector with the cluster indexes of the common cells
% com_clusters = morpho_clusters(ib);
% %and get the data for the corresponding cells
% com_cells = invitro_norm(ia,:);
% 
% %allocate memory for the average profiles
% com_ave = zeros(morpho_clunum,size(invitro_norm,2));
% %allocate memory for the numbers of cells in each cluster
% 
% %for all the morpho_clusters
% for clu = 1:morpho_clunum
%     %skip the morpho_clusters that don't have representatives in the imported
%     %cells
%     if sum(com_clusters==clu)==0
%         continue
%     end
%     %get the common cells for this cluster
%     com_ave(clu,:) = mean(com_cells(com_clusters==clu,:));
% end
% 
% %plot the averages
% figure
% imagesc(com_ave)
% xlabel('PCs')
% ylabel('invitro clusters/morpho numbers')
% set(gca,'XTick',1:size(com_ave,2),'XTickLabels',invitro_vars)
% set(gca,'YTick',1:invitro_clunum,'YTickLabels');
%% Look at the intersection of the clusters, including maps
close all

%get the common cells between in vivo and in vitro
[~,ia,ib] = intersect(invitro_names,morpho_names);

%get the properties of the common cells
invitro_subraw = invitro_raw(ia,:);

%get the morpho subset that is common
morpho_subclusters = morpho_clusters(ib);
%create a color map based on the number of morpho/invitro clusters
morpho_color = jet(morpho_clunum);

%get the invitro cluster indexes of the iviv cells
invitro_subclusters = invitro_clusters(ia);
%get the range for soma plotting
soma_cent = invitro_raw(:,6);
soma_lim = [min(soma_cent)-1,max(soma_cent)+1];
%also get the side bias
side_bias = invitro_raw(:,7);

%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invitro_clunum
    %if the cluster only contains 1 cell
    if sum(invitro_clusters==clu)==1
        %skip it
        continue
    end
    %average the cells in question
    clu_cells = squeeze(mean(invitro_maps(:,:,invitro_clusters==clu,:),3));
    
    %create a figure for the combined plot
    h = figure;
    %plot the overlap map
    subplot(1,2,1)
    map_plot2(clu_cells,strcat('Cluster No:',num2str(clu),',Members:',num2str(sum(invitro_clusters==clu))),1,h,1,1)
    %get the depths in the cluster
    soma_clu = soma_cent(invitro_clusters==clu);
    %extract the side bias
    side_clu = side_bias(invitro_clusters==clu,:);
    %also plot the distribution of depths in the cluster
    subplot(1,2,2)
    %first the responses
    plot(side_clu,soma_clu,'ok')
    hold('on')
    
    %then the average
    errorbar(mean(side_clu(:,1)),mean(soma_clu),std(soma_clu)./sqrt(length(soma_clu)),'go','MarkerFaceColor','g')
    
    %counter of morpho cells
    morpho_count = 0;
    %for all the invivo invitro_clusters
    for fclu = 1:morpho_clunum
        %get the cells from this invivo cluster and plot
        soma_funcclu = invitro_subraw(morpho_subclusters==fclu&invitro_subclusters==clu,6);
        side_funcclu = invitro_subraw(morpho_subclusters==fclu&invitro_subclusters==clu,7);
        
        plot(side_funcclu,soma_funcclu,'*','MarkerFaceColor',morpho_color(fclu,:))
        %count the morpho cells
        morpho_count = morpho_count + length(side_funcclu);
    end
    set(gca,'YLim',soma_lim,'XLim',[-1.1 1.1],'Ydir','reverse')
    %plot a cross in 0,0
    plot(zeros(2,1),get(gca,'YLim'),'-k')
    plot(get(gca,'XLim'),[sum(soma_lim)/2,sum(soma_lim)/2],'-k')
    yyaxis right
    ylabel(strcat('Soma depth,','Morpho cells:',num2str(morpho_count)))
    xlabel('Side bias')
    
    set(gca,'YTick',[],'YColor','k')

end
%% Produce matrix plots of interaction between clusters
close all

%get the common cells between in vivo and in vitro
[~,ia,ib] = intersect(invitro_names,morpho_names);

% %get the properties of the common cells
% invitro_subraw = invitro_raw(ia,:);

%get the morpho subset that is common
morpho_subclusters = morpho_clusters(ib);

%get the invitro cluster indexes of the iviv cells
invitro_subclusters = invitro_clusters(ia);

%allocate memory for the matrix
interaction_mat = zeros(invitro_clunum,morpho_clunum);

%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invitro_clunum
 
    %for all the invivo invitro_clusters
    for fclu = 1:morpho_clunum
        %calculate the overlap between the two groups
        interaction_mat(clu,fclu) = sum(morpho_subclusters==fclu&invitro_subclusters==clu);
    end
end

%plot the resulting matrix
figure
imagesc(interaction_mat)
xlabel('Morpho clusters')
ylabel('Invitro clusters')
set(gca,'XTick',1:morpho_clunum)
set(gca,'YTick',1:invitro_clunum)
colorbar

%get the common cells between in vivo and in vitro
[~,ia,ib] = intersect(invitro_names,invivo_names);

% %get the properties of the common cells
% invitro_subraw = invitro_raw(ia,:);

%get the morpho subset that is common
invivo_subclusters = invivo_clusters(ib);

%get the invitro cluster indexes of the iviv cells
invitro_subclusters = invitro_clusters(ia);

%allocate memory for the matrix
interaction_mat = zeros(invitro_clunum,invivo_clunum);

%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:invitro_clunum
 
    %for all the invivo invitro_clusters
    for fclu = 1:invivo_clunum
        %calculate the overlap between the two groups
        interaction_mat(clu,fclu) = sum(invivo_subclusters==fclu&invitro_subclusters==clu);
    end
end

%plot the resulting matrix
figure
imagesc(interaction_mat)
xlabel('Invivo clusters')
ylabel('Invitro clusters')
set(gca,'XTick',1:invivo_clunum)
set(gca,'YTick',1:invitro_clunum)
colorbar

%get the common cells between in vivo and in vitro
[~,ia,ib] = intersect(morpho_names,invivo_names);

% %get the properties of the common cells
% morpho_subraw = morpho_raw(ia,:);

%get the morpho subset that is common
invivo_subclusters = invivo_clusters(ib);

%get the invitro cluster indexes of the iviv cells
morpho_subclusters = morpho_clusters(ia);

%allocate memory for the matrix
interaction_mat = zeros(morpho_clunum,invivo_clunum);

%plot maps of a set of invitro_clusters
%for all the invitro_clusters
for clu = 1:morpho_clunum
 
    %for all the invivo invitro_clusters
    for fclu = 1:invivo_clunum
        %calculate the overlap between the two groups
        interaction_mat(clu,fclu) = sum(invivo_subclusters==fclu&morpho_subclusters==clu);
    end
end

%plot the resulting matrix
figure
imagesc(interaction_mat)
xlabel('Invivo clusters')
ylabel('Morpho clusters')
set(gca,'XTick',1:invivo_clunum)
set(gca,'YTick',1:morpho_clunum)
colorbar
%% Plot the PCs from all datasets



figure
subplot(1,2,1)
imagesc(invitro_coeff)
set(gca,'YTick',1:size(invitro_norm,2),'YTickLabels',invitro_vars)
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
subplot(1,2,2)
imagesc(invitro_cluave)
set(gca,'XTick',1:size(invitro_norm,2),'XTickLabels',invitro_vars,'XTickLabelRotation',90)
set(gca,'YTick',1:invitro_clunum,'YTickLabels',invitro_clumem)
set(gca,'TickLabelInterpreter','none')
ylabel('Cluster')


figure
subplot(1,2,1)
imagesc(invivo_coeff)
set(gca,'YTick',1:size(invivo_norm,2),'YTickLabels',invivo_vars)
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
subplot(1,2,2)
% %blank the last cluster cause it only has one member and messes up the
% %scale
% invivo_cluave(3,:) = 0;
imagesc(invivo_cluave)
set(gca,'XTick',1:size(invivo_norm,2),'XTickLabels',invivo_vars,'XTickLabelRotation',90)
set(gca,'YTick',1:invivo_clunum,'YTickLabels',invivo_clumem)
set(gca,'TickLabelInterpreter','none')
ylabel('Cluster')

figure
subplot(1,2,1)
imagesc(morpho_coeff)
set(gca,'YTick',1:size(morpho_norm,2),'YTickLabels',morpho_vars)
set(gca,'TickLabelInterpreter','none')
xlabel('PCs')
subplot(1,2,2)
imagesc(morpho_cluave)
set(gca,'XTick',1:size(morpho_norm,2),'XTickLabels',morpho_vars,'XTickLabelRotation',90)
set(gca,'YTick',1:morpho_clunum,'YTickLabels',morpho_clumem)
set(gca,'TickLabelInterpreter','none')
ylabel('Cluster')
%% Plot variables from the clusters

close all

%define the target data set
tar_set = invitro_raw;
tar_vars = invitro_vars;
tar_clunum = invitro_clunum;
tar_cluind = invitro_clusters;
tar_clus = [2 3 4];
%define the target variables
% tar_par = [5 6 3 4 1 12];
tar_par = 1:7;

%get the sub number of clusters
sub_clunum = length(tar_clus);
%define the colors
color_vec = jet(sub_clunum);
%get the number of plots
plot_num = length(tar_par);

figure

plot_matrix = tar_set(:,tar_par);
plot_matrix(plot_matrix(:,1:4)==0) = NaN;
%for all the clusters
for clu = 1:tar_clunum
    if any(tar_clus==clu)
        continue
    end
    tar_cluind(tar_cluind==clu) = NaN;
end

%for all the plots
for plots = 1:plot_num
    subplot(round(sqrt(plot_num)),ceil(sqrt(plot_num)),plots)
    
     boxplot(plot_matrix(:,tar_par(plots)),tar_cluind,'Colors',color_vec,'Notch','on','PlotStyle','traditional')
%     %counter for the x position
%     x_count = 1;
%     %for each selected cluster
%     for clu = tar_clus
%         curr_clu = plot_matrix(tar_cluind==clu,:);        
%         plot((rand(size(curr_clu,1),1)./2)-0.25 + x_count,curr_clu(:,tar_par(plots)),'*')
%         hold('on')
%         errorbar(x_count,nanmean(curr_clu(:,tar_par(plots))),nanstd(curr_clu(:,tar_par(plots)))./sqrt(size(curr_clu,1)),'*')
%         %update the x counter
%         x_count = x_count + 1;
%     end
    ylabel(tar_vars{plots})
    xlabel('Clusters')
    set(gca,'XTick',1:sub_clunum,'XTickLabels',1:sub_clunum)

end
%% Statistical tests

