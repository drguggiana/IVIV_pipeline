str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% Define array names and get the arrays without NaNs etc.
field_list = {'excMap','inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
iviv_cells = [str(:).iviv]==1;
morpho_cells = ~cellfun(@isempty, {str.morph});
% cell_idx = find(iviv_cells&morpho_cells);
% cell_idx = find(iviv_cells);
% cell_idx = find(morpho_cells);
cell_idx = 1:length(str);
% get the number of cells to include
cell_num = length(cell_idx);
% allocate memory for the individual cells
cell_cell = cell(cell_num,1);
% for all the cells
for cells = 1:cell_num
    % allocate a small cell to hold each field
    field_cell = cell(field_number,1);
    field_cell_raw = cell(field_number,1);
    
    % for all the fields
    for fields = 1:field_number
        field_cell{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'excMap')
            field_cell{fields} = field_cell{fields}./min(field_cell{fields});
        elseif contains(field_list{fields}, 'inhMap')
            field_cell{fields} = field_cell{fields}./max(field_cell{fields});
        end       
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});    
    
    for fields = 1:field_number
        field_cell_raw{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'excMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        elseif contains(field_list{fields}, 'inhMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        end       
    end
    
    cell_cell_raw{cells} = vertcat(field_cell_raw{:});
end
% concatenate the results
cell_cell = cat(2,cell_cell{:})';
cell_cell_raw = cat(2,cell_cell_raw{:})';
% remove cells with NaNs
non_nan_cells = sum(isnan(cell_cell),2)==0;
nan_vector = find(non_nan_cells>0);
cell_cell = cell_cell(non_nan_cells,:);
cell_cell_raw = cell_cell_raw(non_nan_cells,:);
% copy the cell to have the original maps later
original_maps = cell_cell;
original_maps_raw = cell_cell_raw;
% redefine cell number based on the rows that didn't contain NaN
cell_num = size(cell_cell,1);
%Pia vector for ex and inh maps 
pia_input=original_maps(:,end);
incl_idx=1;
%% Get 16x16 maps for ex and in
ex_map = reshape(original_maps(:,1:256)',16,16,length(nan_vector));
in_map = reshape(original_maps(:,257:512)',16,16,length(nan_vector));
% Get 16x16 maps for ex and in RAW
ex_map_raw = reshape(original_maps_raw(:,1:256)',16,16,length(nan_vector));
in_map_raw = reshape(original_maps_raw(:,257:512)',16,16,length(nan_vector));
%% %% Calculate overlap of the maps ex and inh
%for all the cells
for cells = 1:length(nan_vector(incl_idx:end));
    %get the binarized maps
    bin_exc = str(nan_vector(cells)).excMap(3:16,:);
    bin_excm= [zeros(2,16); bin_exc];
    bin_inh = str(nan_vector(cells)).inhMap;
    %get the layers
   % lyr = invitro_struct(cells).layers;
    %get the overlap map
    ov_map(:,:,cells) = squeeze(sum(cat(3,bin_excm,bin_inh),3));  
    ov_map_n(:,:,cells)=ov_map(:,:,cells)./max(max(ov_map(:,:,cells)));
    ex_tot(:,cells)=str(nan_vector(cells)).excinhTotal(1);
    in_tot(:,cells)=str(nan_vector(cells)).excinhTotal(2);
    %layers(:,:,cells)=str(nan_vector(i)).layers;
end


%% total sum per map
for i=1:length(nan_vector);
    temp=ov_map_n(:,:,i);
    sum_ex(i)=sum(temp(:));
    temp2=in_map(:,:,i);
    sum_in(i)=sum(temp2(:));
    temp=[];
    temp2=[];
end
%% Morphology and Cell ID and iviv cell ID
%morphtraces are always there but str.morph are the ones that are decent
%traced and used for analysis
for i=1:length(nan_vector)
    cellID_str(i)=str(nan_vector(i)).cellID;
    if ~isempty(str(nan_vector(i)).morph)==1;
        morph_cells(i)=1;
        m_flip_a(i)=str(nan_vector(i)).morph_flip_again;
    else
        morph_cells(i)=0;
        m_flip_a(i)=NaN;
    end
end
morph_cells_id=find(morph_cells==1);
iviv_celid=find([str(nan_vector).iviv]==1);
%% Get the morphology density maps
for i=1:length(nan_vector)
    if m_flip_a(i)==1
        morpho_basal{i,:}=fliplr(str(nan_vector(i)).morphoMap_basal);
        morpho_apical{i,:}=fliplr(str(nan_vector(i)).morphoMap_apical);
    else
       morpho_basal{i,:}=str(nan_vector(i)).morphoMap_basal;
       morpho_apical{i,:}=str(nan_vector(i)).morphoMap_apical;
    end
end

%% Align the exc and inh maps vertically
% distance between squares
grid_distance = 69;
% get the square edges
edges = 0:grid_distance:16*grid_distance;
% get the exc and inh maps
% clean_ex = reshape(cat(3,str(non_nan_cells).excMap),256,[]);
% clean_in = reshape(cat(3,str(non_nan_cells).inhMap),256,[]);
clean_ex = cell_cell(:,1:256)';
clean_in = cell_cell(:,257:512)';
clean_exraw = cell_cell_raw(:,1:256)';
clean_inraw = cell_cell_raw(:,257:512)';
pia_bins = discretize(cell_cell(:,513),edges);
clean_ov=reshape(ov_map,256,148);
clean_diff=reshape(diff_map,256,148);
% get the number of actually populated bins
bin_num = max(pia_bins);
% allocate memory for the aligned maps
aligned_maps_ex = zeros(cell_num,16+bin_num,16);
aligned_maps_in = zeros(cell_num,16+bin_num,16);
aligned_maps_exraw = zeros(cell_num,16+bin_num,16);
aligned_maps_inraw = zeros(cell_num,16+bin_num,16);
aligned_maps_ov = zeros(cell_num,16+bin_num,16);
aligned_maps_diff = zeros(cell_num,16+bin_num,16);
aligned_maps_basal = zeros(cell_num,16+bin_num,16);
aligned_maps_apical = zeros(cell_num,16+bin_num,16);

% for all the cells
for cells = 1:cell_num
    % get the cells displacement
    cell_soma = pia_bins(cells);   
    % place the map into the aligned matrix based on this displacement
    aligned_maps_ex(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_ex(:,cells),16,16);
    aligned_maps_in(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_in(:,cells),16,16);
    aligned_maps_ov(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_ov(:,cells),16,16);
    aligned_maps_diff(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_diff(:,cells),16,16);
    aligned_maps_exraw(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_exraw(:,cells),16,16);
    aligned_maps_inraw(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_inraw(:,cells),16,16);
    % if the morpho is missing, skip
    if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
        aligned_maps_basal(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
            morpho_basal{cells};        
    end
    % if the morpho is missing, skip
    if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
        aligned_maps_apical(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
            morpho_apical{cells};        
    end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_ov = reshape(aligned_maps_ov,cell_num,[]);
aligned_maps_diff = reshape(aligned_maps_diff,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
aligned_maps_basal = reshape(aligned_maps_basal,cell_num,[]);
aligned_maps_apical = reshape(aligned_maps_apical,cell_num,[]);
%% Align the exc and inh maps horizontally
clean_ex = aligned_maps_ex';
clean_in = aligned_maps_in';
clean_ov = aligned_maps_ov';
clean_diff = aligned_maps_diff';
clean_basal = aligned_maps_basal;
clean_apical = aligned_maps_apical;
% get the center of the soma
soma_centers = cat(1,str(:).somaCenter);
soma_centers = soma_centers(non_nan_cells,1)- min(soma_centers(non_nan_cells,1));
% % define the cell number
% cell_num = size(clean_ex,2);
% bin the pial vector
% pia_bins = discretize(cell_cell(:,513),bin_num);
% center_bins = discretize(soma_centers, hbin_num);
center_bins = discretize(soma_centers, edges);
% if all the cells are in the same bin, don't center
if length(unique(center_bins)) == 1
    hbin_num = 0;
else
    % get the number of actually populated bins
    hbin_num = max(center_bins);
    % allocate memory for the aligned maps
    aligned_maps_ex = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_in = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_exraw = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_inraw = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_ov = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_diff = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_basal = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_apical = zeros(cell_num,16+bin_num, 16+hbin_num);
    % for all the cells
    for cells = 1:cell_num
        % get the cells displacement
        cell_soma = center_bins(cells);   
        % place the map into the aligned matrix based on this displacement
        aligned_maps_ex(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_ex(:,cells),16+bin_num,16);
        aligned_maps_in(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_in(:,cells),16+bin_num,16);
        aligned_maps_exraw(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_exraw(:,cells),16+bin_num,16);
        aligned_maps_inraw(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_inraw(:,cells),16+bin_num,16);
        aligned_maps_ov(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_ov(:,cells),16+bin_num,16);
        aligned_maps_diff(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_diff(:,cells),16+bin_num,16);
        % if the morpho is missing, skip
        if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
            aligned_maps_basal(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_basal(cells,:,:);
        end
        % if the morpho is missing, skip
        if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
            aligned_maps_apical(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_apical(cells,:,:);
        end
    end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
aligned_maps_ov = reshape(aligned_maps_ov,cell_num,[]);
aligned_maps_diff = reshape(aligned_maps_diff,cell_num,[]);
aligned_maps_basal = reshape(aligned_maps_basal,cell_num,[]);
aligned_maps_apical = reshape(aligned_maps_apical,cell_num,[]);


%% Run 2 separate PCAs on ex and in and basal and apical morpho
%which maps to include: 
%incl_idx=65;
[coeff_ex,score_ex,latent_ex,~,explained_ex,mu] = pca(aligned_maps_ex(incl_idx:end,:,:));
[coeff_in,score_in,latent_in,~,explained_in,mu] = pca(aligned_maps_in(incl_idx:end,:,:));
[coeff_exraw,score_exraw,latent_exraw,~,explained_exraw,mu] = pca(aligned_maps_exraw(incl_idx:end,:,:));
[coeff_inraw,score_inraw,latent_inraw,~,explained_inraw,mu] = pca(aligned_maps_inraw(incl_idx:end,:,:));
[coeff_ov,score_ov,latent_ov,~,explained_ov,mu] = pca(aligned_maps_ov(incl_idx:end,:,:));
[coeff_diff,score_diff,latent_diff,~,explained_diff,mu] = pca(aligned_maps_diff(incl_idx:end,:,:));
[coeff_com,score_com,latent_com,~,explained_com,mu] = pca([cat(1,aligned_maps_ex) cat(1,aligned_maps_in)]);
[coeff_basal,score_basal,latent_basal,~,explained_basal,mu] = pca(aligned_maps_basal(morph_cells_id,:,:));
[coeff_apical,score_apical,latent_apical,~,explained_apical,mu] = pca(aligned_maps_apical(morph_cells_id,:,:));
pia_input=pia_input(incl_idx:end);
%% Display variance explained for ex and in maps
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y2=explained_in;yyaxis left;h2=bar(y2);h2.EdgeColor = 'k';h2.FaceColor = [0 0 0.7];%ylim([0 20])
hold on;y1=explained_ex;yyaxis left;h=bar(y1);h.EdgeColor = 'k';h.FaceColor = [0.7 0 0];
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:5:25]);
y3=cumsum(y1);y4=cumsum(y2);yyaxis right
%cdfplot(r.explained);hold on;
p1=plot(y3,'--');hold on;p2=plot(y4,'--')
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[1 0 0];p2.Color=[0 0 1];legend('Inhibition', 'Excitation','Excitation', 'Inhibition');%legend box off;
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]);
%xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','Pial depth'})
%% Display variance explained for apical and basal density
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y2=explained_apical;yyaxis left;h2=bar(y2);h2.EdgeColor = 'k';h2.FaceColor = [0 0 0.7];%ylim([0 20])
hold on;y1=explained_basal;yyaxis left;h=bar(y1);h.EdgeColor = 'k';h.FaceColor = [0.7 0 0];
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:5:100]);
y3=cumsum(y1);y4=cumsum(y2);yyaxis right
%cdfplot(r.explained);hold on;
p1=plot(y3,'--');hold on;p2=plot(y4,'--')
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[1 0 0];p2.Color=[0 0 1];legend('Apical', 'Basal','Basal', 'Apical');%legend box off;
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]);
%% %% Display variance explained for overlap
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y1=explained_ov;yyaxis left;h=bar(y1);h.EdgeColor = 'k';h.FaceColor = [0.7 0 0];
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:5:100]);
y3=cumsum(y1);
%cdfplot(r.explained);hold on;
p1=plot(y3,'--');hold on;
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[1 0 0];
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]);

%% Display coefficent of PCs ALIGNED in and ex
fig2 = figure;
set(fig2, 'Name', 'PCs coeff ex in aligned maps');
set(fig2, 'Position', [0, 0, 800, 300]);set(gcf,'color','w')
set(gcf,'color','w')
[cmap]=buildcmap('kwg');
%[cmap]=buildcmap('kwb');
for i = 1:3
    subplot(2,3,i)
    imagesc(reshape(coeff_ex(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.2 0.4]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
hold on;

for i = 1:3
    subplot(2,3,i+3)
    imagesc(reshape(coeff_in(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.2 0.4]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
%Use the scores from first 3 PCs and add togéther 
data_w_input=[score_ex(:,1:3) score_in(:,1:3)];
%% %% Display coefficent of PCs ALIGNED basal and apical
fig2 = figure;
set(fig2, 'Name', 'PCs aligned maps');
set(fig2, 'Position', [0, 0, 800, 300]);set(gcf,'color','w')
set(gcf,'color','w')
[cmap]=buildcmap('kwg');
for i = 1:3
    subplot(2,3,i)
    imagesc(reshape(coeff_apical(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.35 0.65]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
hold on;

for i = 1:3
    subplot(2,3,i+3)
    imagesc(reshape(coeff_basal(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.35 0.65]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
%Use the scores from first 3 PCs and add together 
data_w_morph=[score_apical(:,1:3) score_basal(:,1:3)];
%% %% Display coefficent of PCs ALIGNED overlap/diff
fig2 = figure;
set(fig2, 'Name', 'PCs aligned maps');
set(fig2, 'Position', [0, 0, 800, 300]);set(gcf,'color','w')
set(gcf,'color','w')
[cmap]=buildcmap('kwg');
for i = 1:3
    subplot(1,3,i)
    imagesc(reshape(coeff_diff(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.35 0.65]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
hold on;

c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
%Use the scores from first 3 PCs and add together 
data_w_ov=[score_ov(:,1:3)];
%% %% Display coefficent of PCs ALIGNED COMBINED
fig2 = figure;
set(fig2, 'Name', 'PCs coeff ex in aligned maps');
set(fig2, 'Position', [0, 0, 1000, 500]);set(gcf,'color','w')
for i = 1:3
    subplot(2,3,i)
    imagesc(reshape(coeff_com(442:end,i),16+bin_num,16+hbin_num));
    axis square;colorbar;colormap('hot');
end

%Use the scores from first 3 PCs and add togéther 
data_w_input_com=[score_com(:,1:3)];

%% Histogram of pia distribution with color coded in vivo morpho by itself
figure;set(gcf,'color','w');
histogram(pia_input,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
hold on;histogram(pia_input(morph_cells_id),'FaceColor','m','Orientation','horizontal');
hold on;histogram(pia_input(iviv_celid),'FaceColor','g','Orientation','horizontal');
legend([' Input (n=' num2str(length(pia_input)),')'],[' Morph. (n=' num2str(length(pia_input(morph_cells_id))),')']...
    ,[' In vivo (n=' num2str(length(pia_input(iviv_celid))),')']);
legend boxoff  
set(gca,'Ydir','reverse')
%histogram(pia_input,'FaceColor','k','FaceAlpha',0.1)
set(gca,'FontSize',14);
%set(gca,'FontWeight','bold')
%% Display correlation between all PCs and pial depth Multiple comparison 
com=[data_w_input(morph_cells_id,:) data_w_morph pia_input(morph_cells_id)]; 
correlation_matrix(com,1);
hold on;
xticks([1:1:13]);yticks([1:1:13]);
xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','PC1_{ap}','PC2_{ap}','PC3_{ap}','PC1_{ba}','PC2_{ba}','PC3_{ba}','Pial depth'});
% names={'PC1ex','PC2ex','PC3ex','PC1in','PC2in','PC3in','Pial depth'};
% set(gca,'XTickLabel',names,'FontSize',10); 
xtickangle(45);
% set(gca,'YTickLabel',names,'FontSize',10); 
ytickangle(45);
yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','PC1_{ap}','PC2_{ap}','PC3_{ap}','PC1_{ba}','PC2_{ba}','PC3_{ba}','Pial depth'})
%% Strongest correlation with pial depth, which is Pc1ex and PC2in
figure;set(gcf,'color','w');
 [R P]=corrcoef(data_w_input(:,1),pia_input);scatter(data_w_input(:,1),pia_input,25,'o','MarkerEdgeColor','k');box off;xlabel('PC1ex');ylabel('Pial depth (µm)');
 text(1.8,120,['r= ' mat2str(round(R(2),2))]);xlim([-2 2]);set(gca,'Ydir','reverse');
 if P(2)<0.05 & P(2)>0.01
     text(1.8,140,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
     text(1.8,140,['p<0.01'])
 elseif P(2)<0.001
     text(1.8,140,['p<0.001']);
 else
     text(1.8,140,['n.s']);
 end
  P = polyfit(data_w_input(:,1),pia_input,1);
    yfit = P(1)*data_w_input(:,1)+P(2);
    %figure;plot(data_w_input(:,1),pia_input,'ko');
    hold on;box off;xlabel('PC1ex','Color','r');ylabel('Pial depth (µm)');
    plot(data_w_input(:,1),yfit,'g-');set(gca,'box','off');set(gcf,'color','w');axis square;
 
 figure;set(gcf,'color','w');
 [R P]=corrcoef(data_w_input(:,5),pia_input);scatter(data_w_input(:,5),pia_input,25,'o','.','MarkerEdgeColor','k');box off;xlabel('PC1ex');ylabel('Pial depth (µm)');
 text(1.8,120,['r= ' mat2str(round(R(2),2))]);xlim([-2 2]);set(gca,'Ydir','reverse');
 if P(2)<0.05 & P(2)>0.01
     text(1.8,140,['p<0.05'])

 elseif P(2)<0.01 & P(2)>0.001
      text(1.8,140,['p<0.01'])
 elseif P(2)<0.001
     text(1.8,140,['p<0.001']);
     else
     text(1.8,140,['n.s']);
 end
  P = polyfit(data_w_input(:,5),pia_input,1);
    yfit = P(1)*data_w_input(:,5)+P(2);
    %figure;plot(data_w_input(:,1),pia_input,'ko');
    hold on;box off;xlabel('PC2in','Color','b');ylabel('Pial depth (µm)');set(gca,'Ydir','reverse');
    plot(data_w_input(:,5),yfit,'g-');set(gca,'box','off');set(gcf,'color','w');axis square;
%% Display input data
%Load fraction and absolute
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv pialD layer_assign] = iviv_profiles(nan_vector(incl_idx:end),str);

%% Calculate the fraction from the overlap for verti and hori
[frac_ovh abs_ovh frac_ovv abs_ovv] = iviv_profiles_ov(ov_map);

%% Calculate the maximum hozontal span
for i=1:length(nan_vector(incl_idx:end))
tmp=find(frac_exh(i,:)>0);
ex_spanh(i)=tmp(end)-tmp(1);
tmp=[];
tmp=find(frac_inh(i,:)>0);
in_spanh(i)=tmp(end)-tmp(1);
tmp=[];
end
%% 
for i=1:length(nan_vector(incl_idx:end))
diffL23(i)=sum(sum(diff_map(3:5,:,i),2));
diffL4(i)=sum(sum(diff_map(6:8,:,i),2));
diffL5(i)=sum(sum(diff_map(9:12,:,i),2));
end
%% %%Display
frac_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv];
abs_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv];
fig7= figure;
set(fig7, 'Name', 'Input distribution');
set(fig7, 'Position', [200, 0, 400, 520]);
set(gcf,'color','w');
subplot(2,2,1);
hold on;
for i=1:length(nan_vector(incl_idx:end))
exp=plot(frac_exv_m(i,:)'*-1,1:16,'-r');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(nanmean(frac_exv_m)'*-1,1:16,nanstd(frac_exv_m)'/sqrt(length(nan_vector)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
%xticklabels({'1','3','5','7','9','11','13','15'});
hold on;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'})
%Inhbition vert
hold on;
for i=1:length(nan_vector(incl_idx:end))
exp=plot(frac_inv(i,:)',1:16,'-b');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(nanmean(frac_inv)',1:16,nanstd(frac_inv)'/sqrt(length(nan_vector)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;
%xlim([0 1]);
hold on;
ylim([1 16]);yticks([1:5:16])
yticklabels({'1','6','11','16'});
%overlap verti
subplot(2,2,3);
hold on;
% for i=1:length(nan_vector(incl_idx:end))
% exp=plot(frac_ovv(i,:)',1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
% exp.Color(4) = 0.2;
% end
% hold on;
% mexp=errorbar(nanmean(frac_ovv)',1:16,nanstd(frac_ovv)'/sqrt(length(nan_vector)),'horizontal','k');set(gca,'Ydir','reverse');
tmp=frac_exv_m(i,:)'-frac_inv(i,:)';
for i=1:length(nan_vector(incl_idx:end))
exp=plot(tmp,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(tmp,1:16,tmp/sqrt(length(nan_vector)),'horizontal','k');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-0.6 0.6]);xticks([-0.6:0.3:0.6]);xticklabels({'-0.6','-0.3','0','0.3','0.6'});
hold on;
ylim([1 16]);yticks([1:5:16])
yticklabels({'1','6','11','16'})
%Excitation hori
subplot(2,2,2);
hold on;
for i=1:length(nan_vector(incl_idx:end))
exp=plot(frac_exh(i,:)'*-1,'-r');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(nanmean(frac_exh)*-1,nanstd(frac_exh)/sqrt(length(nan_vector)),'k');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
%Inhbition hori
hold on;
for i=1:length(nan_vector(incl_idx:end))
exp=plot(frac_inh(i,:)','-b');;xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.2;
exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(nanmean(frac_inh),nanstd(frac_inh)/sqrt(length(nan_vector)),'k');
mexp.CapSize=3;
%ylim([0 1]);
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'})
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');
%overlap hori
subplot(2,2,4);
hold on;
for i=1:length(nan_vector(incl_idx:end))
exp=plot(frac_ovh(i,:)','-m');;xlabel('Horizontal input');ylabel('Fraction of total input');box off;

exp.Color(4) = 0.2;
end
hold on;
mexp=errorbar(nanmean(frac_ovh),nanstd(frac_ovh)/sqrt(length(nan_vector)),'k');
mexp.CapSize=3;
ylim([-0.6 0.6]);yticks([-0.6:0.3:0.6]);yticklabels({'-0.6','-0.3','0','0.3','0.6'});
hold on;
xlim([1 16]);;xticks([1:5:16])
xticklabels({'1','6','11','16'});hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');

%% Get 16x16 maps for basal and apical
for i=1:length(morph_cells_id)
     if isempty(morpho_basal{morph_cells_id(i),:,:})==0
     ba_map(:,:,i)=morpho_basal{morph_cells_id(i),:,:};
     ap_map(:,:,i)=morpho_apical{morph_cells_id(i),:,:};
     else
     ba_map(:,:,i)=ones(16,16)*NaN;
     ap_map(:,:,i)=ones(16,16)*NaN;
     end
 end
%% Plot average maps for all cells
F = figure;
set(gcf,'color','w');
set(F, 'Name', 'Correlation pial depth');
set(F, 'Position', [200, 0, 800, 200]);
subplot(1,3,1);
sf=1;    
exc_map =nanmean(ex_map,3);
inh_map =nanmean(in_map,3);
ove_map = cat(3,exc_map,inh_map);
%define the plot type (2 for excitatory)
explot_type = 2;
%get the pial distance
soma_info = pia_input(1);
map_plot3(exc_map,'',explot_type,F,sf,0,1);
hold on;        
      x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
    %plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3);
   
    hold on;
  end
  subplot(1,3,2);
  inplot_type = 3;
 map_plot3(inh_map,'',inplot_type,F,sf,0,1);
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
   % plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3)
    hold on;
  end
  subplot(1,3,3);
   bplot_type = 1;
map_plot3(ove_map,'',bplot_type,F,sf,0,1)
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
   % plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3)
    hold on;
  end
  %% Plot individual maps (choose any)
  %plot_maps(str,ex_map,in_map,nan_vector,10,pia_input);
  plot_maps(str,ex_map_raw,in_map_raw,nan_vector,95,pia_input);
   % cells with ori pref > 30: 65    80    85    88    91    93    98   100   101   104   112   119   122   126   131   138   142   148
 %   cells with ori pref < 30: 66    68    70    72    78    79    84    86    87    95    97   102   103   105   106   107   108   113   114   115   120   121   123   127   128   129   130   133   135   136   139   140   141
 % 143   144   145

%% Clustering HCA VERTICAL
clu_num = 3;
%pcs =[];
pcs     =[1 2 3 4 5 6];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%% Ex and in avergae maps
F = figure;
set(F, 'Name', 'Heat average input maps for cluster');
set(F, 'Position', [0, 0, 1000, 1000]);set(gcf,'color','w')
 set(gcf,'color','w');
sf=1;
for clu = 1:clu_num
    p_i=[0 nanmean(pia_input(find(idx_input==clu)))];   
       excc_map =nanmean(ex_map(:,:,idx_input==clu),3)
       inhc_map =nanmean(in_map(:,:,idx_input==clu),3)
       ove_map = cat(3,excc_map,inhc_map);
         %define the plot type (2 for excitatory)
        explot_type = 2;
         %get the pial distance
        soma_info = pia_input(idx_input==clu);
        subplot(3,3,clu)
        map_plot3(excc_map,'',explot_type,F,sf,1,1);
             title(num2str(sum(idx_input==clu)));
        inplot_type = 3;
        subplot(3,3,clu+3);
        map_plot3(inhc_map,'',inplot_type,F,sf,1,1);
  
   bplot_type = 1;
   subplot(3,3,clu+6);
   map_plot3(ove_map,'',bplot_type,F,sf,1,1);
    p_i=[];
    subplot(3,3,clu+3)
end
sf=1;
%% Display parameters for clusters: INPUT Horizontal
%data_input=[frac_exv(:,1:9) frac_inv(:,1:9)];
data_input=[];
data_input=[frac_exh(incl_idx:end,4:13) frac_inh(incl_idx:end,4:13)];
input_clusterh=cluster_show(idx_input,[data_input pia_input],clu_num,2);
%% Display parameters for clusters: INPUT Vertical
data_input=[];
data_input=[frac_exv(incl_idx:end,3:9) frac_inv(incl_idx:end,1:9)];
%data_input=[frac_exh(:,4:13) frac_inh(:,4:13)];
input_clusterv=cluster_show(idx_input,[data_input pia_input],clu_num,2);
%% %% Display parameters for clusters: PCs
data_input=[];
data_input=[data_w_input(incl_idx:end,:)];
%data_input=[frac_exh(:,4:13) frac_inh(:,4:13)];
input_clusterpc=cluster_show(idx_input,[data_input pia_input],clu_num,2);
%% %% %%Display
frac_exv_c1=frac_exv_m(find(idx_input==1),:);
frac_exv_c2=frac_exv_m(find(idx_input==2),:);
frac_exv_c3=frac_exv_m(find(idx_input==3),:);
fig20= figure;
set(fig20, 'Name', 'Input distribution');
set(fig20, 'Position', [200, 0, 1000, 1040]);
set(gcf,'color','w');
subplot(2,2,1);
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*-1,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c2)
exp=plot(frac_exv_c2(i,:)'*-1,1:16,'-y');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*-1,1:16,'-g');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_exv_c1)'*-1,1:16,nanstd(frac_exv_c1)'/sqrt(length(frac_exv_c1)),'horizontal','m');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
hold on;
mexp=errorbar(nanmean(frac_exv_c2)'*-1,1:16,nanstd(frac_exv_c2)'/sqrt(length(frac_exv_c2)),'horizontal','y');set(gca,'Ydir','reverse');
mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_exv_c3)'*-1,1:16,nanstd(frac_exv_c3)'/sqrt(length(frac_exv_c3)),'horizontal','g');set(gca,'Ydir','reverse');
mexp.CapSize=3;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'});

frac_exv_c1=frac_inv(find(idx_input==1),:);
frac_exv_c2=frac_inv(find(idx_input==2),:);
frac_exv_c3=frac_inv(find(idx_input==3),:);
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*1,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c2)
exp=plot(frac_exv_c2(i,:)'*1,1:16,'-y');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*1,1:16,'-g');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_exv_c1)'*1,1:16,nanstd(frac_exv_c1)'/sqrt(length(frac_exv_c1)),'horizontal','m');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
hold on;
mexp=errorbar(nanmean(frac_exv_c2)'*1,1:16,nanstd(frac_exv_c2)'/sqrt(length(frac_exv_c2)),'horizontal','y');set(gca,'Ydir','reverse');
mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_exv_c3)'*1,1:16,nanstd(frac_exv_c3)'/sqrt(length(frac_exv_c3)),'horizontal','g');set(gca,'Ydir','reverse');
mexp.CapSize=3;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'});

frac_exv_c1=frac_ovv(find(idx_input==1),:);
frac_exv_c2=frac_ovv(find(idx_input==2),:);
frac_exv_c3=frac_ovv(find(idx_input==3),:);
%overlap verti
subplot(2,2,3);
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*1,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c2)
exp=plot(frac_exv_c2(i,:)'*1,1:16,'-y');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*1,1:16,'-g');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_exv_c1)'*1,1:16,nanstd(frac_exv_c1)'/sqrt(length(frac_exv_c1)),'horizontal','m');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
hold on;
mexp=errorbar(nanmean(frac_exv_c2)'*1,1:16,nanstd(frac_exv_c2)'/sqrt(length(frac_exv_c2)),'horizontal','y');set(gca,'Ydir','reverse');
mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_exv_c3)'*1,1:16,nanstd(frac_exv_c3)'/sqrt(length(frac_exv_c3)),'horizontal','g');set(gca,'Ydir','reverse');
mexp.CapSize=3;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'});%


frac_exv_c1=frac_exh(find(idx_input==1),:);
frac_exv_c2=frac_exh(find(idx_input==2),:);
frac_exv_c3=frac_exh(find(idx_input==3),:);
subplot(2,2,2);
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*-1,'-m');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c1)*-1,nanstd(frac_exv_c1)/sqrt(length(frac_exv_c1)),'m');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c2(i,:)'*-1,'-y');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c2)*-1,nanstd(frac_exv_c2)/sqrt(length(frac_exv_c2)),'y');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*-1,'-g');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c3)*-1,nanstd(frac_exv_c3)/sqrt(length(frac_exv_c3)),'g');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');


frac_exv_c1=frac_inh(find(idx_input==1),:);
frac_exv_c2=frac_inh(find(idx_input==2),:);
frac_exv_c3=frac_inh(find(idx_input==3),:);
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*1,'-m');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c1)*1,nanstd(frac_exv_c1)/sqrt(length(frac_exv_c1)),'m');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c2(i,:)'*1,'-y');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c2)*1,nanstd(frac_exv_c2)/sqrt(length(frac_exv_c2)),'y');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*1,'-g');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c3)*1,nanstd(frac_exv_c3)/sqrt(length(frac_exv_c3)),'g');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--')

frac_exv_c1=frac_ovh(find(idx_input==1),:);
frac_exv_c2=frac_ovh(find(idx_input==2),:);
frac_exv_c3=frac_ovh(find(idx_input==3),:);
subplot(2,2,4);
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c1(i,:)'*1,'-m');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c1)*1,nanstd(frac_exv_c1)/sqrt(length(frac_exv_c1)),'m');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c1)
exp=plot(frac_exv_c2(i,:)'*1,'-y');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c2)*1,nanstd(frac_exv_c2)/sqrt(length(frac_exv_c2)),'y');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;
for i=1:length(frac_exv_c3)
exp=plot(frac_exv_c3(i,:)'*1,'-g');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_c3)*1,nanstd(frac_exv_c3)/sqrt(length(frac_exv_c3)),'g');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;line([8 8], [-1 1],'Color','k','LineStyle','--')

%% Difference in pial depth of clusters
pia_cluster=pia_input;
comi=[pia_cluster(find(idx_input==1)); pia_cluster(find(idx_input==2));  pia_cluster(find(idx_input==3))]';
group=[ones(1,length(pia_cluster(find(idx_input==1)))) ones(1,length(pia_cluster(find(idx_input==2))))*2 ones(1,length(pia_cluster(find(idx_input==3))))*3];
figure; set(gcf,'color','w');
plot(group(1:length(pia_cluster(find(idx_input==1)))),comi(1:length(pia_cluster(find(idx_input==1))))','mo');
%hold on;plot(group(length(pia_cluster(find(idx_input==2)))+1:end),comi(length(oric(find(cluster_id==3)))+1:end)','mo');
hold on;plot(ones(1,length(pia_cluster(find(idx_input==2))))*2 ,pia_cluster(find(idx_input==2)),'yo');
hold on;plot(ones(1,length(pia_cluster(find(idx_input==3))))*3 ,pia_cluster(find(idx_input==3)),'go');
%plotSpread(comi'-90,'categoryIdx',group, 'categoryMarkers',{'.','.'},'categoryColors',{'k','m'});
hold on;
boxplot(comi,group,'Colors','k','OutlierSize',0.001);
box off;xticklabels({'Cluster 1','Cluster 2','Cluster 3'});ylabel('Pial depth');
set(gca,'Ydir','reverse');
%% Assemble the feature vector
 cell_cell = cat(2,score_ex(:,1:3),score_in(:,1:3));
 %% Extract Oripref, Dirpref and sigmatuning for ipsi, contra and bino cell. For bino cell I use the orientation preference for the 
%more dominant eye (based on ODI)
non_nan_idx=nan_vector(1:end);
for i=1:length(non_nan_idx);
    idx_amp_inf(i,:)=str(non_nan_idx(i)).hori_peak_pl;
if ~isempty(str(non_nan_idx(i)).Ori)==1 & isempty(str(non_nan_idx(i)).unres)==1
    %ori_both(:,i)=str(non_nan_idx(i)).ori_a;
    iv_ODI(i)=str(non_nan_idx(i)).ODI;
    if str(non_nan_idx(i)).contra==1
     oripref(i)=str(non_nan_idx(i)).Ori(1); 
     dirpref(i)=str(non_nan_idx(i)).Dir(1);
     sigma(i)=str(non_nan_idx(i)).sigma(1);
     fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(1:360);
     iv_OSI(i)=str(non_nan_idx(i)).OSI(1);
     iv_DSI(i)=str(non_nan_idx(i)).DSI(1);
     iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(1);
    elseif str(non_nan_idx(i)).ipsi==1;
      oripref(i)=str(non_nan_idx(i)).Ori(2);
      dirpref(i)=str(non_nan_idx(i)).Dir(2);
      sigma(i)=str(non_nan_idx(i)).sigma(2);
      fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(361:end);
      iv_OSI(i)=str(non_nan_idx(i)).OSI(2);
      iv_DSI(i)=str(non_nan_idx(i)).DSI(2);
      iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(2);
    elseif str(non_nan_idx(i)).bino==1;
         if str(non_nan_idx(i)).ODI>=0;
        oripref(i)=str(non_nan_idx(i)).Ori(1);
        dirpref(i)=str(non_nan_idx(i)).Dir(1);
        sigma(i)=str(non_nan_idx(i)).sigma(1);
        fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(1:360);
        iv_OSI(i)=str(non_nan_idx(i)).OSI(1);
        iv_DSI(i)=str(non_nan_idx(i)).DSI(1);
        iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(1);
         else str(non_nan_idx(i)).ODI<0;
            oripref(i)=str(non_nan_idx(i)).Ori(2);
        dirpref(i)=str(non_nan_idx(i)).Dir(2);
        sigma(i)=str(non_nan_idx(i)).sigma(2);
        fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(361:end);
        iv_OSI(i)=str(non_nan_idx(i)).OSI(2);
        iv_DSI(i)=str(non_nan_idx(i)).DSI(2);
        iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(2);
    end
 else str(non_nan_idx(i)).unres==1;
        oripref(i)=NaN;
        dirpref(i)=NaN;
        sigma(i)=NaN;
        fit_resp(i,:)=ones(360,1)*NaN;
        iv_ODI(i)=NaN;
        iv_OSI(i)=NaN;
        iv_DSI(i)=NaN;
        iv_Ca(i)=NaN;
        %ori_both(:,i)=ones(1,2)*NaN;
 end

else
oripref(i)=NaN;
dirpref(i)=NaN;
sigma(i)=NaN;
fit_resp(i,:)=ones(360,1)*NaN;
iv_ODI(i)=NaN;
iv_OSI(i)=NaN;
iv_DSI(i)=NaN;
iv_Ca(i)=NaN;
end
end
%% Extract other in vivo paramters: PCI and spontenaeous events 
for i=1:length(non_nan_idx);   
if ~isempty(str(non_nan_idx(i)).pci)==1
iv_spon(i)=str(non_nan_idx(i)).sad;
iv_popcop(i)=str(non_nan_idx(i)).pci;
else
iv_spon(i)=NaN;    
iv_popcop(i)=NaN; 
end
end
%% %OD decompistion
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).iv_OD_decom)==1
        od_decomp(i,:)=str(non_nan_idx(i)).iv_OD_decom;
    else
        od_decomp(i,:)=ones(1,24)*NaN;
    end
end
% SF decomposition
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).iv_SFTF_decom)==1
        sftf_decomp(i,:)=str(non_nan_idx(i)).iv_SFTF_decom;
    else
        sftf_decomp(i,:)=ones(1,24)*NaN;
    end
end
%% SFTF both eyes
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).sftf_resp)==1
        iv_SF(i)=str(non_nan_idx(i)).SF;
        iv_TF(i)=str(non_nan_idx(i)).TF;
        oripref_sftf(i)=str(non_nan_idx(i)).Ori_sftf;
        dirpref_sftf(i)=str(non_nan_idx(i)).Dir_sftf;
        oripref_sftf_all(i,:)=str(non_nan_idx(i)).Ori_sftf_all;
        dirpref_sftf_all(i,:)=str(non_nan_idx(i)).Dir_sftf_all;
        osi_sftf(i)=str(non_nan_idx(i)).OSI_sftf;
        dsi_sftf(i)=str(non_nan_idx(i)).DSI_sftf;
        osi_sftf_all(i,:)=str(non_nan_idx(i)).OSI_sftf_all;
        dsi_sftf_all(i,:)=str(non_nan_idx(i)).DSI_sftf_all;
        if sum(str(non_nan_idx(i)).sftf_resp)==0;
            iv_resp2(i)=0;
        else sum(str(non_nan_idx(i)).sftf_resp)>0;
            iv_resp2(i)=1;
        end
    else
        iv_SF(i)=NaN;
        iv_TF(i)=NaN;  
        iv_resp2(i)=NaN;
        oripref_sftf(i)=NaN;
        dirpref_sftf(i)=NaN;
        osi_sftf(i)=NaN;
        dsi_sftf(i)=NaN;
        osi_sftf_all(i,:)=ones(1,9)*NaN;
        dsi_sftf_all(i,:)=ones(1,9)*NaN;
        oripref_sftf_all(i,:)=ones(1,9)*NaN;
        dirpref_sftf_all(i,:)=ones(1,9)*NaN;
    end
end
%% 
%Extract hemisphere and ori
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).sliceOri)==1
        slice_ori(i)=str(non_nan_idx(i)).sliceOri;
       slice_hemi(i)=str(non_nan_idx(i)).hemisphere;
    else
        slice_ori(i)=NaN;
        slice_hemi(i)=NaN;      
    end
end
% Extract morpho parameters

for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).morph)==1
        morph_matrix(i,:)=str(non_nan_idx(i)).morph;
      
    else
        morph_matrix(i,:)=ones(1,24)*NaN;
      
    end
end
% Extract ipsi contra bino
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).ipsi)==1
        iv_ipsi(i,:)=str(non_nan_idx(i)).ipsi;
      
    else
        iv_ipsi(i,:)=NaN;
      
    end
end
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).contra)==1
        iv_contra(i,:)=str(non_nan_idx(i)).contra;
      
    else
        iv_contra(i,:)=NaN;
      
    end
end
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).bino)==1
        iv_bino(i,:)=str(non_nan_idx(i)).bino;
      
    else
        iv_bino(i,:)=NaN;
      
    end
end
%% Clusters with in vivo paramaters
%which clusters out 
clun=2;
cluster_id=idx_input;
cluster_id(find(idx_input==clun))=[];
oric=oripref-90;
%oric=abs(od_decomp(:,2)')
pia=pia_input;
%oric=pia';
pia(find(idx_input==clun))=[];
oric(find(idx_input==clun))=[];
%% Clusters bar plot to show differences 
comi=[oric(find(cluster_id==1)) oric(find(cluster_id==3))];
group=[ones(1,length(oric(find(cluster_id==1)))) ones(1,length(oric(find(cluster_id==3))))*2];
figure; set(gcf,'color','w');
plot(group(1:length(oric(find(cluster_id==1)))),comi(1:length(oric(find(cluster_id==1))))','mo');
hold on;plot(group(length(oric(find(cluster_id==3)))+1:end),comi(length(oric(find(cluster_id==3)))+1:end)','go');
%plotSpread(comi'-90,'categoryIdx',group, 'categoryMarkers',{'.','.'},'categoryColors',{'k','m'});
hold on;
boxplot(comi,group,'Colors','k','OutlierSize',0.001);
box off;xticklabels({'Cluster 1','Cluster 3'});ylabel('Orientation preference');
c1_ori=oripref(find(idx_input==1));
c3_ori=oripref(find(idx_input==3));
ylim([-90 90]);
yticks(-90:22.5:90)
%yticklabels({'-90','6','11','16'});

%% %% Display parameters for clusters: in vivo
data_input=[];
data_input=[oripref'  dirpref' iv_OSI' iv_DSI' iv_ODI' iv_SF' iv_TF' iv_popcop' iv_spon'];
%data_input=[frac_exh(:,4:13) frac_inh(:,4:13)];
input_clusteriviv=cluster_show(idx_input,[data_input pia_input],clu_num,2);

%% Simple Correlation matrix between input, invivo , morphology
%only use iviv cells
scores_test=[diffL23; diffL4; diffL5]';
idxiv=find(isnan(oripref)==0);
com_iviv=[scores_test(idxiv,:) pia_input(idxiv,:) ex_spanh(idxiv)' in_spanh(idxiv)' iv_spon(idxiv)' iv_popcop(idxiv)' iv_DSI(idxiv)' iv_OSI(idxiv)' iv_ODI(idxiv)' oripref(idxiv)' iv_Ca(idxiv)' sigma(idxiv)']; 
%com_iviv=[data_w_input(idxiv,:) pia_input(idxiv,:) iv_spon(idxiv)' iv_popcop(idxiv)' iv_DSI(idxiv)' iv_OSI(idxiv)' iv_ODI(idxiv)' oripref(idxiv)' dirpref(idxiv)' iv_Ca(idxiv)', oripref_sftf(idxiv)']; 
correlation_matrix(zscore(com_iviv),0);
hold on;
xticks([1:1:16]);
% xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','Pial depth','SAD (Hz)','PCI','gDSI','gOSI','ODI','Ori (°)','Dir (°)','pCa','SF','TF','ori2'});
% % names={'PC1ex','PC2ex','PC3ex','PC1in','PC2in','PC3in','Pial depth'};
% % set(gca,'XTickLabel',names,'FontSize',10); 
% xtickangle(45);
% % set(gca,'YTickLabel',names,'FontSize',10); 
% ytickangle(45);
% yticks([1:1:16]);
% yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','Pial depth','SAD (Hz)','PCI','gDSI','gOSI','ODI','Ori (°)','Dir (°)','pCa','SF','TF'});
%% %% Simple Correlation matrix between input and SF TF
scores_test=[diffL23; diffL4; diffL5]';
idxiv=[];
idxiv=find(isnan(oripref_sftf)==0);
com_iviv=[scores_test(idxiv,:) oripref_sftf(idxiv)' dirpref_sftf(idxiv)' iv_TF(idxiv)' iv_SF(idxiv)' osi_sftf(idxiv)' dsi_sftf(idxiv)']; 
%com_iviv=[data_w_input(idxiv,:) pia_input(idxiv,:) iv_spon(idxiv)' iv_popcop(idxiv)' iv_DSI(idxiv)' iv_OSI(idxiv)' iv_ODI(idxiv)' oripref(idxiv)' dirpref(idxiv)' iv_Ca(idxiv)', oripref_sftf(idxiv)']; 
correlation_matrix(zscore(com_iviv),0);
hold on;
xticks([1:1:16]);
%% %% Simple Correlation matrix between input and SF TF
scores_test=[score_ex(:,1:3) score_in(:,1:3)] ;
idxiv=[];
idxiv=find(isnan(oripref_sftf_all(:,7))==0);
com_iviv=[scores_test(idxiv,:)  oripref_sftf_all(idxiv,7)]; 
%com_iviv=[data_w_input(idxiv,:) pia_input(idxiv,:) iv_spon(idxiv)' iv_popcop(idxiv)' iv_DSI(idxiv)' iv_OSI(idxiv)' iv_ODI(idxiv)' oripref(idxiv)' dirpref(idxiv)' iv_Ca(idxiv)', oripref_sftf(idxiv)']; 
correlation_matrix(zscore(com_iviv),0);
hold on;
xticks([1:1:16]);
%% Morphology correlation 
idxmo=find(isnan(morph_matrix(:,1))==0);
idxiv=find(isnan(oripref)==0);
[a b c]=intersect(idxmo,idxiv);
%com_iviv=[data_w_input(idxiv,:) pia_input(idxiv,:),morph_matrix(idxiv,:)]; 
% com_iviv=[morph_matrix(idxiv,:) iv_spon(idxiv)' iv_popcop(idxiv)' iv_DSI(idxiv)' iv_OSI(idxiv)' iv_ODI(idxiv)' oripref(idxiv)' dirpref(idxiv)' iv_peakCa(idxiv)' iv_SF(idxiv)' iv_TF(idxiv)']; 
com_iviv=[morph_matrix(a,:) iv_spon(a)' iv_popcop(a)' iv_DSI(a)' iv_OSI(a)' iv_ODI(a)' oripref(a)' dirpref(a)' iv_Ca(a)']; 
correlation_matrix(zscore(com_iviv),0);
hold on;
%% 
a=[];
a=morph_cells_id;
com_iviv=[data_w_morph(:,:) iv_spon(a)' iv_popcop(a)' iv_DSI(a)' iv_OSI(a)' iv_ODI(a)' oripref(a)' dirpref(a)' iv_Ca(a)']; 
correlation_matrix(com_iviv,0);
%% Correlations PC scores and in vivo parameters, plot the significant coorelation between PC2in and Ori
par1=cell_cell(:,5);
%par2=oripref-90;
par2=oripref_sftf-90;
scoren='PC2_{in}';
figure;scatter(par1,par2,'ko');set(gcf,'color','w');
idxn = find(~isnan(par2));
P = polyfit(par1(idxn),par2(idxn)',1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:m
ax(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','b');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(1.1,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(1.1,70,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
      text(1.1,70,['p<0.01'])
 elseif P(2)<0.001
     text(1.1,70,['p<0.001'])
 else
     text(1.1,70,['n.s.'])
 end
 xlim([-2 2]);
%% Correlations PC scores and in vivo parameters, plot the significant coorelation between PC2in and Ori
par1=cell_cell(:,5);
%par2=dpref-90;
par2=oripref-90;
par3=nan_vector;
%par3=idx_input;
idf=1:length(pia_input);
scoren='PC2_{in}';
[cmap]=buildcmap('kmgy');
pointsize=20;
figure;scatter(par1(idf),par2(idf),pointsize,par3(idf),'filled');
 set(gcf,'color','w');
 box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (\mum)'; c.Ticks=[min(par3):round(min(par3)/3,-1):max(par3)];
 set( c, 'YDir', 'reverse' );
 caxis([min(par3) max(par3)]);
idxn = find(~isnan(par2));
P = polyfit(par1(idxn),par2(idxn)',1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','b');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(1.1,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(1.1,70,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
      text(1.1,70,['p<0.01'])
 elseif P(2)<0.001
     text(1.1,70,['p<0.001'])
 else
     text(1.1,70,['n.s.'])
 end
 xlim([-2 2]);
 %axis square
 %% to check cell IDS
 figure;scatter3(par1(idf),par2(idf),par3(idf),'filled');
 %% %% Correlations PC scores and in vivo parameters + CLUSTERS
 par1=cell_cell(:,5);
%par2=dpref-90;
par2=oripref-90;
%par3=pia_input;
par3=idx_input;
idf=1:length(pia_input);
scoren='PC2_{in}';
[cmap]=buildcmap('gm');
pointsize=20;
figure;scatter(par1(idf),par2(idf),pointsize,par3(idf),'filled');
 set(gcf,'color','w');
 %box off; c=colorbar;
 colormap(cmap);
%  c.Label.String = 'Cluster'; c.Ticks=[min(par3):round(min(par3)/3,-1):max(par3)];
%  set( c, 'YDir', 'reverse' );
%  caxis([min(par3) max(par3)]);
idxn = find(~isnan(par2));
P = polyfit(par1(idxn),par2(idxn)',1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','b');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(1.1,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(1.1,70,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
      text(1.1,70,['p<0.01'])
 elseif P(2)<0.001
     text(1.1,70,['p<0.001'])
 else
     text(1.1,70,['n.s.'])
 end
 xlim([-2 2]);
 %% Calculate the peak inhbition horizontally
[peak_h jjt]=find(frac_exh'==max(frac_exh'));
[peak_hi jjt]=find(frac_inh'==max(frac_inh'));
%% Get the aligned maps
for i=1:length(nan_vector)
align_ex(:,:,i)=reshape(aligned_maps_ex(i,:),16+bin_num,16+hbin_num);
align_in(:,:,i)=reshape(aligned_maps_in(i,:),16+bin_num,16+hbin_num);
end
%% Display input data
[frac_inha abs_inha frac_ovvaa abs_ovvaa] = iviv_profiles_ov(align_in);
[frac_exha abs_exha frac_ovvaa abs_ovvaa] = iviv_profiles_ov(align_ex);

[peak_hia jjt]=find(frac_inha'==max(frac_inha'));
[peak_ha jjt]=find(frac_exha'==max(frac_exha'));

%%  horizontal peak for ex and in  per row
%EX
for i=1:length(frac_exh)
temp_map=ex_map(3:16,:,i);
counter=1;
for k=1:14
frac_hpv(k,i,:)=temp_map(k,:)./sum(temp_map);
counter=counter+1;
end
end
%IN
for i=1:length(frac_exh)
temp_map=in_map(1:16,:,i);
counter=1;
for k=1:16
frac_hpvi(k,i,:)=temp_map(k,:)./sum(temp_map);
counter=counter+1;
end
end
%%  horizontal peak for ex and in  per row
mape=align_ex;
mapi=align_in;
for i=1:length(frac_exh)
temp_map=mape(:,:,i);
counter=1;
for k=3:14
idx_map=find(temp_map(k,:)./sum(temp_map(k,:))==max(temp_map(k,:)./sum(temp_map(k,:))));
if ~isempty(idx_map)==1
idx_map_ex(i,counter)=idx_map(1);
else
    idx_map_ex(i,counter)=NaN;
end
counter=counter+1;
end

end
for i=1:length(frac_exh)
temp_map=mapi(:,:,i);
counter=1;
for k=3:14
idx_map=find(temp_map(k,:)./sum(temp_map(k,:))==max(temp_map(k,:)./sum(temp_map(k,:))));
if ~isempty(idx_map)==1
idx_map_in(i,counter)=idx_map(1);
else
    idx_map_in(i,counter)=NaN;
end
counter=counter+1;
end

end
%peak_diff_lay=idx_map_ex-idx_map_in;
%% Get indices of visual responsive cells amnong the 148 cells
id_ori=find(~isnan(oripref));

for i=1:length(nan_vector)
somax(i)=str(nan_vector(i)).somaCenter(1);
somay(i)=552-str(nan_vector(i)).somaCenter(2);
end
%% Calculate centroid of each map: FOR ENTIRE MAP
%call function centroid_map
%IN
[sx sy wx_in wy_in ang_ai ang_bi ang_ai_v hypo_i vec_slope_i qdi] = centroid_map(in_map(:,:,:),somax,pia_input,[1:148],0);
[sx sy wx_inl4 wy_inl4 ang_ail4 ang_bil4 ang_ai_vl4 hypo_il4 vec_slope_il4 qdil4] = centroid_map(in_map(6:8,:,:),somax,pia_input,[1:148],5);
[sx sy wx_inl3 wy_inl3 ang_ail3 ang_bil3 ang_ai_vl3 hypo_il3 vec_slope_il3 qdil3] = centroid_map(in_map_raw(3:5,:,:),somax,pia_input,[1:148],2);
%EX
[sx sy wx_ex wy_ex ang_ae ang_be ang_ae_v hypo_e vec_slope_e qde] = centroid_map(ex_map,somax, pia_input,[1:148],0);
[sx sy wx_exl4 wy_exl4 ang_ael4 ang_bel4 ang_ae_vl4 hypo_el4 vec_slope_el4 qdel4] = centroid_map(ex_map(6:8,:,:),somax,pia_input,[1:148],5);
[sx sy wx_exl3 wy_exl3 ang_ael3 ang_bel3 ang_ae_vl3 hypo_el3 vec_slope_el3 qdel3] = centroid_map(ex_map(3:5,:,:),somax,pia_input,[1:148],2);
%OV
[sx sy wx_ov wy_ov ang_ov ang_bov ang_ov_v hypo_ov vec_slope_ov qdov] = centroid_map(ov_map,somax, pia_input,[1:148],0);
[sx sy wx_ovl4 wy_ovl4 ang_ovl4 ang_ovbl4 ang_ov_vl4 hypo_ovl4 vec_slope_ovl4 qdovl4] = centroid_map(ov_map(6:8,:,:),somax,pia_input,[1:148],5);
[sx sy wx_ovl3 wy_ovl3 ang_ovl3 ang_ovbl3 ang_ov_vl3 hypo_ovl3 vec_slope_ovl3 qdovl3] = centroid_map(ov_map(3:5,:,:),somax,pia_input,[1:148],2);
%APICAL
[sx_m sy_m wx_ap wy_ap ang_aap ang_bap ang_aap_v hypo_ap vec_slope_ap qda] = centroid_map(ap_map,somax,pia_input,morph_cells_id,0);
%BASAL
[sx_m sy_m wx_ba wy_ba ang_aba ang_bba ang_aba_v hypo_ba vec_slope_ba gdb] = centroid_map(ba_map,somax,pia_input,morph_cells_id,0);
%% Plot vectors for all INHIBITION
figure;hold on;set(gcf,'color','w');
dpx=[];dpy=[];
% if 8-wx<0
%    wxtemp=
% else
% end
dpx=wx_inl3-sx;
dpy=wy_inl3-sy;
id_p=[1:148];
for t=1:length(id_p)
q=quiver(8,sy(t),dpx(t),dpy(t),0);
%q=quiver(c_cord(t,1),4,dpy(t),dpx(t),0);
 q.Color='cyan';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
hold on;plot(8,sy(t),'c^');
%hold on;plot(c_cord(t,1),4,'k^');
end
hold on;
for t=1:length(id_ori);
q=quiver(8,sy(t),dpx(id_ori(t)),dpy(id_ori(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(8,sy(t),'k^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='black';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([1 16]);ylim([1 16]);

%% Plot vectors for all Excitation 
figure;hold on;set(gcf,'color','w');
dpx=[];dpy=[];
dpx=wx_ex-sx;
dpy=wy_ex-sy;
id_p=[1:148];
for t=1:length(id_p)
q=quiver(8,3,dpx(t),dpy(t),0);
%q=quiver(c_cord(t,1),4,dpy(t),dpx(t),0);
 q.Color='red';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
hold on;plot(8,3,'r^');
%hold on;plot(c_cord(t,1),4,'k^');
end
hold on;
for t=1:length(id_ori);
q=quiver(8,3,dpx(id_ori(t)),dpy(id_ori(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(8,3,'k^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='black';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([1 16]);ylim([1 16]);
%% Plot vectors for all apical and basal
figure;hold on;set(gcf,'color','w');
dpx=[];dpy=[];
dpx=wx_ap-sx_m;
dpy=wy_ap-sy_m;
id_p=[1:98];
for t=1:length(id_p)
q=quiver(sx(t),sy(t),dpx(t),dpy(t),0);
%q=quiver(c_cord(t,1),4,dpy(t),dpx(t),0);
 q.Color='red';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
hold on;plot(sx(t),sy(t),'b^');
%hold on;plot(c_cord(t,1),4,'k^');
end
 hold on;
 dpx=[];dpy=[];
dpx=wx_ba-sx_m;
dpy=wy_ba-sy_m;
id_p=[1:98];
for t=1:length(id_p)
q=quiver(sx(t),sy(t),dpx(t),dpy(t),0);
%q=quiver(c_cord(t,1),4,dpy(t),dpx(t),0);
 q.Color='black';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
hold on;plot(sx(t),sy(t),'b^');
%hold on;plot(c_cord(t,1),4,'k^');
end
% for t=1:length(id_ori);
% q=quiver(sx(id_ori(t)),sy(id_ori(t)),dpy(id_ori(t)),dpx(id_ori(t)),0);
% %q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
% hold on;plot(sx(id_ori(t)),sy(id_ori(t)),'k^');
% %plot(c_cord(id_ori(t),1),4,'g^');
% q.Color='black';
%  q.MaxHeadSize=0.1;
%  q.LineWidth=1;
% end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([4 15]);ylim([1 10]);
%% %% Plot vectors for all INHIBITION
m=4;
color_v=[1 0 0;0 1 0; 0 0 1; 0 0 0];
[dis_ori E]=discretize((oripref(id_ori)-90),m);

dpx=[];dpy=[];
dpx=wx_inl3-sx;
dpy=wy_inl3-sy;
id_p=[1:148];
idxm1=find(dis_ori==1);
idxm2=find(dis_ori==2);
idxm3=find(dis_ori==3);
idxm4=find(dis_ori==4);
figure;polarplot(oripref(id_ori(idxm1)),ang_ail3(id_ori(idxm1)),'ro')
hold on;
polarplot(oripref(id_ori(idxm2)),ang_ail3(id_ori(idxm2)),'go')
hold on;
polarplot(oripref(id_ori(idxm3)),ang_ail3(id_ori(idxm3)),'bo')
hold on;
polarplot(oripref(id_ori(idxm4)),ang_ail3(id_ori(idxm4)),'ko')
figure;hold on;
for k=1:m
idxm=find(dis_ori==k);
hold on;
for t=1:length(idxm);
q=quiver(8,3,dpx(id_ori(idxm(t))),dpy(id_ori(idxm(t))),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(8,3,'k^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color=color_v(k,:);
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
idxm=[];
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([1 16]);ylim([1 16]);
%% single example with vector
t=10;
inh_map=in_map(:,:,id_ori(t));
F = figure;
set(gcf,'color','w');
inplot_type = 3;
sf=1;
map_plot3(inh_map,'',inplot_type,F,sf,0,1);
hold on;plot(sx(id_ori(t)),sy(id_ori(t)),'k^')
hold on;set(gcf,'color','w');q=quiver(sx(id_ori(t)),sy(id_ori(t)),(wx_in(id_ori(t))-sx(id_ori(t))),(wy_in(id_ori(t))-sy(id_ori(t))),0);
 q.Color='black';
 q.MaxHeadSize=0.8;
 q.LineWidth=1;
 title([num2str(ang_ai(id_ori(t))) ',' num2str(data_w_input(id_ori(t),5))]);
 
% xlim([1 16]);ylim([1 16]);
% set(gca,'Ydir','reverse');
%% %% Correlations PC scores and in vivo parameters, plot the significant coorelation between PC2in and Ori
par1=data_w_input(id_ori,5);
par2=gg';
scoren='PC2_{in}';
figure;scatter(par1,par2,'ko');set(gcf,'color','w');
idxn = find(~isnan(par2));
P = polyfit(par1(idxn),par2(idxn)',1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel('Angle towards centre of mass(°)');xlabel(scoren,'Color','b');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(1.8,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(1.8,70,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
      text(1.8,70,['p<0.01'])
 elseif P(2)<0.001
     text(1.8,70,['p<0.001'])
 else
     text(1.8,70,['n.s.'])
 end
 xlim([-2 2]);
 %% 
 avg_ang=nanmean([ang_ail4;ang_ail3])
 %% par
 figure;scatter3(par1,par2,par3)
 %% 
 pia_idx=find(pia_input>200)
 %% 
 %cellv=1:148;
%par1=ang_ail3(id_ori)
par1=diffL23(id_ori);
par2=iv_spon(id_ori);
par3=osi_sftf_all(id_ori,7);
% par1=ang_ail3(id_ori)-ang_ael3(id_ori);
% par2=oripref(id_ori)-90;
% par3=iv_OSI(id_ori);
%hypo_il3 vec_slope_il4
pointsize=20;
scoren='Angle towards centre of mass(°)';
figure;scatter(par1,par2,pointsize,par3,'filled');set(gcf,'color','w');colorbar;
idxn=find(~isnan(par1));
P = polyfit(par1(idxn),par2(idxn),1);
 yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');;p2.Color(4) = 0.5;
 ylabel('Orientation preference (°)');xlabel(scoren,'Color','k');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(90,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(90,70,['p<0.05'])
 elseif P(2)<0.01 & P(2)>0.001
      text(90,70,['p<0.01'])
 elseif P(2)<0.001
     text(90,70,['p<0.001'])
 else
     text(90,70,['n.s.'])
 end
 %xlim([-2 2]);

 %% Outliers
 
 outli=[65 88 104 105 120 122 140 142 148];
 
 figure;hold on
 for i=1:length(outli)
     subplot(3,3,i)
     hold on;
     imagesc(in_map(:,:,outli(i)));colorbar;
     set(gca,'Ydir','reverse');
    
     hold on;plot(sx(outli(i)),sy(outli(i)),'w^');
     hold on;plot(wx_inl3(outli(i)),wy_inl3(outli(i)),'ko');
 end
 %% 
 
 %% Plot correlation between line length and angles

 figure;set(gcf,'color','w');scatter(hypo_e,hypo_i);xlabel('Excitation');ylabel('Inhibition');hold on;refline(1,0);
 figure;set(gcf,'color','w');scatter(ang_ae,ang_ai);xlabel('Excitation');ylabel('Inhibition');hold on;refline(1,0);
 ang_ei=ang_ae-ang_ai;
 hypo_ei=hypo_e-hypo_i;

%% responsive vs non responsive
gf=find(iv_resp2==0);
rf=find(iv_resp2==1);
c_resp=rf(find(iv_Ca(find(iv_resp2==1))>200))
% c_sil=gf(find(iv_resp(find(iv_resp2==0))==0));
c_sil=gf;
figure;hold on;set(gcf,'color','w');
dpx=[];dpy=[];
dpx=wx_inl3-sx;
dpy=wy_inl3-sy;

for t=1:length(id_ori);
q=quiver(sx(id_ori(t)),3,dpx(id_ori(t)),dpy(id_ori(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx(id_ori(t)),3,'g^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='green';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end

for t=1:length(c_sil);
q=quiver(sx(c_sil(t)),3,dpx(c_sil(t)),dpy(c_sil(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx(c_sil(t)),3,'r^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='red';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
xlim([4 13])
%% Plot responsive vs unresponsive iviv
% frac_exv_uresp=[];
% frac_exv_resp=[];
% frac_exv_uresp=frac_exv_m(c_sil,:);
% frac_exv_resp=frac_exv_m(c_resp,:);
 frac_exv_uresp=[];
 frac_exv_resp=[];
 frac_exv_uresp=frac_ovv(c_sil,:);
 frac_exv_resp=frac_ovv(c_resp,:);
fig20= figure;
set(fig20, 'Name', 'Input distribution RS vs NVS');
set(fig20, 'Position', [200, 0, 1000, 1040]);
set(gcf,'color','w');
subplot(2,2,1);
hold on;
for i=1:size(frac_exv_resp,1)
exp=plot(frac_exv_resp(i,:)'*-1,1:16,'-g');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:size(frac_exv_uresp,1)
exp=plot(frac_exv_uresp(i,:)'*-1,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end

mexp=errorbar(nanmean(frac_exv_resp)'*-1,1:16,nanstd(frac_exv_resp)'/sqrt(length(frac_exv_resp)),'horizontal','g');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
hold on;
mexp=errorbar(nanmean(frac_exv_uresp)'*-1,1:16,nanstd(frac_exv_uresp)'/sqrt(length(frac_exv_uresp)),'horizontal','m');set(gca,'Ydir','reverse');
mexp.CapSize=3;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'});

hold on;
frac_exv_uresp=[];
frac_exv_resp=[];
frac_exv_uresp=frac_inv(c_sil,:);
frac_exv_resp=frac_inv(c_resp,:);

hold on;
for i=1:size(frac_exv_resp,1)
exp=plot(frac_exv_resp(i,:)'*1,1:16,'-g');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
for i=1:size(frac_exv_uresp,1)
exp=plot(frac_exv_uresp(i,:)'*1,1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end

mexp=errorbar(nanmean(frac_exv_resp)'*1,1:16,nanstd(frac_exv_resp)'/sqrt(length(frac_exv_resp)),'horizontal','g');set(gca,'Ydir','reverse');
mexp.CapSize=3;
xlim([-1 1]);
xticks([-1:0.5:1])
hold on;
mexp=errorbar(nanmean(frac_exv_uresp)'*1,1:16,nanstd(frac_exv_uresp)'/sqrt(length(frac_exv_uresp)),'horizontal','m');set(gca,'Ydir','reverse');
mexp.CapSize=3;
hold on;
ylim([1 16]);
yticks([1:5:16])
yticklabels({'1','6','11','16'});

% frac_exv_uresp=[];
% frac_exv_resp=[];
% frac_exv_uresp=frac_exh(c_sil,:);
% frac_exv_resp=frac_exh(c_resp,:);
frac_exv_uresp=[];
frac_exv_resp=[];
frac_exv_uresp=frac_ovh(c_sil,:);
frac_exv_resp=frac_ovh(c_resp,:);

subplot(2,2,2);
hold on;
for i=1:size(frac_exv_resp,1)
exp=plot(frac_exv_resp(i,:)'*-1,'-g');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_resp)*-1,nanstd(frac_exv_resp)/sqrt(length(frac_exv_resp)),'g');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});

hold on;
for i=1:size(frac_exv_uresp,1)
exp=plot(frac_exv_uresp(i,:)'*-1,'-m');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_uresp)*-1,nanstd(frac_exv_uresp)/sqrt(length(frac_exv_uresp)),'m');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;

hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');

frac_exv_uresp=[];
frac_exv_resp=[];
frac_exv_uresp=frac_inh(c_sil,:);
frac_exv_resp=frac_inh(c_resp,:);
hold on;
for i=1:size(frac_exv_resp,1)
exp=plot(frac_exv_resp(i,:)'*1,'-g');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_resp)*1,nanstd(frac_exv_resp)/sqrt(length(frac_exv_resp)),'g');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});

hold on;
for i=1:size(frac_exv_uresp,1)
exp=plot(frac_exv_uresp(i,:)'*1,'-m');xlabel('Horizontal input');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
hold on;
mexp=errorbar(nanmean(frac_exv_uresp)*1,nanstd(frac_exv_uresp)/sqrt(length(frac_exv_uresp)),'m');
mexp.CapSize=3;
ylim([-1 1]);yticks([-1:0.5:1])
hold on;
xlim([1 16]);xticks([1:5:16])
xticklabels({'1','6','11','16'});
hold on;

hold on;line([8 8], [-1 1],'Color','k','LineStyle','--');

%% Clusters
% gf=find(iv_resp2==0);
% % c_sil=gf(find(iv_resp(find(iv_resp2==0))==0));
% c_sil=gf;
idx_input(id_ori);
c1_idx=find(idx_input(id_ori)==1);
c3_idx=find(idx_input(id_ori)==3);
sx_sub=sx(id_ori);
sy_sub=sy(id_ori);
figure;hold on;set(gcf,'color','w');
dpx=[];dpy=[];
dpx=wx_in(id_ori)-sx_sub;
dpy=wy_in(id_ori)-sy_sub;

for t=1:length(c1_idx);
q=quiver(sx_sub(c1_idx(t)),3,dpx(c1_idx(t)),dpy(c1_idx(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx_sub(c1_idx(t)),3,'m^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='green';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
hold on;

for t=1:length(c3_idx);
q=quiver(sx_sub(c3_idx(t)),3,dpx(c3_idx(t)),dpy(c3_idx(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx_sub(c3_idx(t)),3,'g^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='magenta';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
%xlim([4 13]);
xlim([5 12]);ylim([1 8]);
%% Clusters Morphology
% gf=find(iv_resp2==0);
% % c_sil=gf(find(iv_resp(find(iv_resp2==0))==0));
% c_sil=gf;
% c1_idx=find(idx_input(id_ori)==1);
% c3_idx=find(idx_input(id_ori)==3);
% sx_sub=sx(id_ori);
% sy_sub=sy(id_ori);
% figure;hold on;set(gcf,'color','w');
% dpx=[];dpy=[];
% dpx=wx_exl3(id_ori)-sx_sub;
% dpy=wy_exl3(id_ori)-sy_sub;
c1_idx=[];
c1_idx=find(idx_input(morph_cells_id)==1);
c2_idx=[];
c2_idx=find(idx_input(morph_cells_id)==2);
c3_idx=[];
c3_idx=find(idx_input(morph_cells_id)==3);
dpx=[];dpy=[];
dpx=wx_ba-sx_m;
dpy=wy_ba-sy_m;
figure;hold on;set(gcf,'color','w');
for t=1:length(c1_idx);
q=quiver(sx_m(c1_idx(t)),3,dpx(c1_idx(t)),dpy(c1_idx(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx_m(c1_idx(t)),3,'g^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='green';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
hold on;

for t=1:length(c3_idx);
q=quiver(sx_m(c3_idx(t)),3,dpx(c3_idx(t)),dpy(c3_idx(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx_m(c3_idx(t)),3,'m^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='magenta';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end

for t=1:length(c2_idx);
q=quiver(sx_m(c2_idx(t)),3,dpx(c2_idx(t)),dpy(c2_idx(t)),0);
%q=quiver(c_cord(id_ori(t),1),4,dpy(id_ori(t)),dpx(id_ori(t)),0);
hold on;plot(sx_m(c2_idx(t)),3,'k^');
%plot(c_cord(id_ori(t),1),4,'g^');
q.Color='black';
 q.MaxHeadSize=0.1;
 q.LineWidth=1;
end
set(gca,'Ydir','reverse');
%set(gca,'Xdir','reverse');
%xlim([4 13]);
%% 


%str_id=find([str(:).iviv]==1);
id_m=nan_vector(morph_cells_id);
idx_input_sub=idx_input(morph_cells_id);
c1=id_m(find(idx_input_sub==1));
c2=id_m(find(idx_input_sub==2));
c3=id_m(find(idx_input_sub==3));
figure;
for i=1:length(c1);
if ~isempty(str(c1(i)).morphtraces)==1
       m=plot_tree(str(c1(i)).morphtraces{1,1},[1 0 0],[0 str(c1(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(c1(i)).morphtraces{1,2},[0 0 0],[0 str(c1(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(c1(i)).morphtraces{1,3},[0 0 1],[0 str(c1(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
else
    m=plot(str(c1(i)).pialD,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
end
 set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
end

figure;
for i=1:length(c2);
if ~isempty(str(c2(i)).morphtraces)==1
       m=plot_tree(str(c2(i)).morphtraces{1,1},[1 0 0],[0 str(c2(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(c2(i)).morphtraces{1,2},[0 0 0],[0 str(c2(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(c2(i)).morphtraces{1,3},[0 0 1],[0 str(c2(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
else
    m=plot(str(c2(i)).pialD,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
end
 set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
end

figure;
for i=1:length(c3);
if ~isempty(str(c3(i)).morphtraces)==1
       m=plot_tree(str(c3(i)).morphtraces{1,1},[1 0 0],[0 str(c3(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(c3(i)).morphtraces{1,2},[0 0 0],[0 str(c3(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(c3(i)).morphtraces{1,3},[0 0 1],[0 str(c3(i)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
else
    m=plot(str(c3(i)).pialD,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
end
 set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
end
%% extract Morphological feature from strcuture
for i=1:length(id_m)
    morph_features(i,:)=str(id_m(i)).morph;
end
%% Morphological features comparison of clusters

input_clu_morpho_feat=cluster_show(idx_input_sub,[morph_features pia_input(morph_cells_id)],3,0);

%% 
par1=ang_ail3(id_ori);
%par2=dpref-90;
par2=oripref(id_ori)-90;
%par2=ang_ail4(id_ori);
%par2=ang_alp(id_ori)
par3=idx_input(id_ori);
%par3=idx_input;
idf=1:length(pia_input);
scoren='PC2_{in}';
[cmap]=buildcmap('kmgy');
pointsize=20;
%figure;scatter(par1(idf),par2(idf),pointsize,par3(idf),'filled');
figure;scatter(par1,par2,pointsize,par3,'filled');
 set(gcf,'color','w');
 box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (\mum)'; c.Ticks=[min(par3):round(min(par3)/3,-1):max(par3)];
 colormap(cmap);
 ylabel('Ori Pref (°)');xlabel('Angle Vector (°)');

%% Morpho density corr basal
% calculate the correlation between morpho density and maps
morphoMaps = morpho_basal;
excMaps = {str(non_nan_cells).excMap};
inhMaps = {str(non_nan_cells).inhMap};
non_nan_names = {str(non_nan_cells).cellName};
correlation_values = cell(cell_num,3);
% for all the cells
for cells = 1:cell_num
    if isempty(morphoMaps{cells})
        correlation_values{cells,3} = 'nanCell';
        continue
    end
    % save the name of the cell
    correlation_values{cells,3} = non_nan_names{cells};
    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = excMaps{cells}*-1;
            case 2
                map = inhMaps{cells};
        end
        % calculate the correlation and store
        %correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
        correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
    end  
end
%% Morpho density corr apical
% calculate the correlation between morpho density and maps
morphoMaps = morpho_apical;
excMaps = {str(non_nan_cells).excMap};
inhMaps = {str(non_nan_cells).inhMap};
non_nan_names = {str(non_nan_cells).cellName};
correlation_values_ap = cell(cell_num,3);
% for all the cells
for cells = 1:cell_num
    if isempty(morphoMaps{cells})
        correlation_values_ap{cells,3} = 'nanCell';
        continue
    end
    % save the name of the cell
    correlation_values_ap{cells,3} = non_nan_names{cells};
    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = excMaps{cells}*-1;
            case 2
                map = inhMaps{cells};
        end
        % calculate the correlation and store
        %correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
        correlation_values_ap{cells,maps} = corr(map(:), morphoMaps{cells}(:));
    end  
end
%% Plot histograms for 
figure;histogram([correlation_values{morph_cells_id,1}]);hold on;histogram([correlation_values_ap{morph_cells_id,1}]);
figure;histogram([correlation_values{morph_cells_id,2}],10);hold on;histogram([correlation_values_ap{morph_cells_id,2}],10);
figure;histogram([correlation_values{morph_cells_id,1}],10);hold on;histogram([correlation_values{morph_cells_id,2}],10);
figure;histogram([correlation_values_ap{morph_cells_id,1}],10);hold on;histogram([correlation_values_ap{morph_cells_id,2}],10);
%% 
fig20= figure;
set(fig20, 'Name', 'MorphINPUT');
set(fig20, 'Position', [200, 0, 800, 800]);
set(gcf,'color','w');
subplot(2,2,1);
cd1=cdfplot([correlation_values{morph_cells_id,1}]);hold on;cd2=cdfplot([correlation_values_ap{morph_cells_id,1}]);box off;grid off;title('Excitation'); 
cd1.Color='c';cd2.Color='k';ylabel('Cumulative');xlabel('Correlation morph and map');legend('basal','apical')
subplot(2,2,2);
cd1=cdfplot([correlation_values{morph_cells_id,2}]);hold on;cd2=cdfplot([correlation_values_ap{morph_cells_id,2}]);box off;grid off;title('Inhibition'); 
cd1.Color='c';cd2.Color='k';ylabel('Cumulative');xlabel('Correlation morph and map');legend('basal','apical')
legend('basal','apical')
subplot(2,2,3);
cd1=cdfplot([correlation_values{morph_cells_id,1}]);hold on;cd2=cdfplot([correlation_values{morph_cells_id,2}]);box off;grid off;title('Basal'); 
cd1.Color='r';cd2.Color='b';ylabel('Cumulative');xlabel('Correlation morph and map');legend('Excitation','Inhibition')
subplot(2,2,4);
cd1=cdfplot([correlation_values_ap{morph_cells_id,1}]);hold on;cd2=cdfplot([correlation_values_ap{morph_cells_id,2}]);box off;grid off;title('Apical'); 
cd1.Color='r';cd2.Color='b';ylabel('Cumulative');xlabel('Correlation morph and map');legend('Excitation','Inhibition')
%% 


figure;set(gcf,'color','w');
[cmap]=buildcmap('kwg');
imagesc(reshape(coeff_in(:,2),16+bin_num,16+hbin_num));
axis square;caxis([-0.2 0.4]);axis off
c=colorbar;

c.Label.String = 'Coefficient';
colormap(cmap);
% q.Color='red';
% q.MaxHeadSize=20;
% q.LineWidth=1;



%[cmap]=buildcmap('kwb');
for i = 1:3
    subplot(2,3,i)
    imagesc(reshape(coeff_ex(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.2 0.4]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
hold on;

for i = 1:3
    subplot(2,3,i+3)
    imagesc(reshape(coeff_in(:,i),16+bin_num,16+hbin_num));
    axis square;caxis([-0.2 0.4]);axis off
end
c=colorbar;colormap(cmap);axis off;
c.Label.String = 'Coefficient';
%% Look at input maps and PCIN2
id_ori=find(~isnan(oripref));
dp=i_cord-c_cord;
sf=1;
F=figure;
set(F, 'Name', 'Maps');
set(F, 'Position', [200, 0, 1500, 1000]);set(gcf,'color','w');

inplot_type = 3;
for t=1:length(id_ori);
subplot(7,8,t)
map_plot3(in_map(:,:,id_ori(t)),'',inplot_type,F,sf,0,0);
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');

hold on;%q=quiver(c_cord(id_ori(t),1),c_cord(id_ori(t),2),dp(id_ori(t),1),dp(id_ori(t),2),0);
% q.Color='red';
% q.MaxHeadSize=20;
% q.LineWidth=1;
%scatter(out(id_ori(t),3),out(id_ori(t),2),'.','MarkerFaceColor','r','MarkerEdgeColor','k');
hold on;plot(c_cord(id_ori(t),1),c_cord(id_ori(t),2),'^')
title([num2str(data_w_input(id_ori(t),5)) ',' num2str(oripref(id_ori(t))-90) ',' num2str(pia_input(id_ori(t)))])
end
%% 
figure;hold on;
for t=1:length(id_ori);
imagesc(test_map(:,:,id_ori(t)));
hold on;scatter(out(id_ori(t),3),out(id_ori(t),2),'o','MarkerFaceColor','r','MarkerEdgeColor','k');
end

%% %% Calculate centroid of each map: for specific layer/row
test_map=in_map;
%test_map=align_in;
for i=1:length(test_map)
    for k=1:16;
props = regionprops(true(size(test_map(k,:,i))), test_map(k,:,i), 'WeightedCentroid');
x_weight_r(i,k,:)=props.WeightedCentroid(1);
y_weight_r(i,k,:)=props.WeightedCentroid(2);
%1st alternative
 cell_cord=[552 pia_input(i)];
 input_cord=[x_weight_r(i,k,:)*69 k*69];
 centr_vector_r=polyfit([cell_cord(1),input_cord(1)],[cell_cord(2),input_cord(2)],1);
 centr_slope_r(i,k,:)=centr_vector_r(1);
    end
end
%% 
for i=1:length(test_map)
l4_cord=[nanmean(x_weight_r(i,7,:)) nanmean(y_weight_r(i,7,:))];
l3_cord=[nanmean(x_weight_r(i,2,:)) nanmean(y_weight_r(i,2,:))];
pass_vec=polyfit([l4_cord(1),l3_cord(1)],[l4_cord(2),l3_cord(2)],1);
pass_slope_r(i)=pass_vec(1);
pass_dis(i)=sqrt((l3_cord(1)-l4_cord(1))^2+(l3_cord(2)-l4_cord(2))^2);
end
%% 
for i=1:length(test_map)
diffll(i)=x_weight(i)-y_weight(i)
end
%% 
for i=1:length(test_map)
map_sum(i)=sum(sum(in_map(:,:,i)));
end

for i=1:length(test_map)
map_nonz(i)=length(nonzeros(in_map(:,:,i)));
end
%% 
for i=1:length(test_map)
temp=find(frac_inh(i,:)>0)
hori_span(i)=temp(end)-temp(1);
end
%% 
id_ori=find(~isnan(oripref));

%% Display Centroid with relation to ORI pref
FJ = figure;
set(FJ, 'Name', 'Per row Horizontal Centroid');
set(FJ, 'Position', [0, 0, 1000, 1000]);set(gcf,'color','w')
hold on;
for i=1:12;
par1=x_weight_r(:,i);
%par2=dpref-90;
par2=oripref-90;
par3=pia_input;
scoren=['ML Centroid_{overlap} R' num2str(i)];
[cmap]=buildcmap('mwg');
pointsize=20;
subplot(3,4,i);
scatter(par1,par2,pointsize,par3,'filled');
xlim([3 12]);
 set(gcf,'color','w');
  box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (\mum)'; c.Ticks=[min(par3):round(min(par3)/2,-1):max(par3)];
idxn = find(~isnan(par2));
idxn2 = find(~isnan(par1));
[aa bb]=intersect(idxn,idxn2);
d1=par1(idxn);
d2=par2(idxn);
P = polyfit(d1(bb),d2(bb)',1);
    yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');axis square;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','k');
 
 [R P]=corrcoef(par1,par2,'row','complete');
 text(3.1,80,['r= ' mat2str(round(R(2),2))],'FontSize', 8);

 if P(2)<0.05 & P(2)>0.01
     text(3.1,60,['p<0.05'],'FontSize', 8);

 elseif P(2)<0.01 & P(2)>0.001
      text(3.1,60,['p<0.01'],'FontSize', 8)
 elseif P(2)<0.001
     text(3.1,60,['p<0.001'],'FontSize', 8)
 else
     text(3.1,60,['n.s.'],'FontSize', 8)
 end
end
%% 
par1=x_weight(:,:);
%par2=dpref-90;
par2=oripref-90;
par3=pia_input;
scoren=['ML Centroid_{overlap} overall'];
[cmap]=buildcmap('mwg');
pointsize=20;
figure;scatter(par1,par2,pointsize,par3,'filled');
%xlim([6.5 10.5]);
 set(gcf,'color','w');
  box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (\mum)'; c.Ticks=[min(par3):round(min(par3)/3,-1):max(par3)];
idxn = find(~isnan(par2));
idxn2 = find(~isnan(par1));
[aa bb]=intersect(idxn,idxn2)
d1=par1(idxn);
d2=par2(idxn);
P = polyfit(d1(bb),d2(bb),1);
    yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');axis square;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','k');
 
 [R P]=corrcoef(par1,par2,'row','complete');
 text(6.7,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(6.7,70,['p<0.05'])

 elseif P(2)<0.01 & P(2)>0.001
      text(6.7,70,['p<0.01'])
 elseif P(2)<0.001
     text(6.7,70,['p<0.001'])
 else
     text(6.7,70,['n.s.'])
 end

%% Fraction horizontal inhbition per row
%% Display Centroid with relation to ORI pref
m=12;
FK = figure;
set(FK, 'Name', 'Per row Horizontal Fraction');
set(FK, 'Position', [0, 0, 1000, 1000]);set(gcf,'color','w')
hold on;
for i=1:12;
par1=frac_hpvi(i,:,m);
%par2=dpref-90;
par2=oripref-90;
par3=pia_input;
scoren=['ML Fraction_{in} R' num2str(i)];
[cmap]=buildcmap('mwg');
pointsize=20;
subplot(3,4,i);
scatter(par1,par2,pointsize,par3,'filled');
xlim([0 1]);
 set(gcf,'color','w');
  box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (\mum)'; c.Ticks=[min(par3):round(min(par3)/2,-1):max(par3)];
idxn = find(~isnan(par2));
idxn2 = find(~isnan(par1));
[aa bb]=intersect(idxn,idxn2);
d1=par1(idxn);
d2=par2(idxn);
P = polyfit(d1(bb),d2(bb),1);
    yfit = P(1)*[min(par1):0.1:max(par1)]+P(2);
 hold on;box off;
 p2=plot([min(par1):0.1:max(par1)],yfit,'--','Color', [0.1 0.1 0.1]);set(gca,'box','off');set(gcf,'color','w');axis square;p2.Color(4) = 0.5;
 ylabel('Ori pref (°)');xlabel(scoren,'Color','r');
 
 [R P]=corrcoef(par1,par2,'row','complete');
 text(0.7,80,['r= ' mat2str(round(R(2),2))],'FontSize', 8);

 if P(2)<0.05 & P(2)>0.01
     text(0.7,60,['p<0.05'],'FontSize', 8);

 elseif P(2)<0.01 & P(2)>0.001
      text(0.7,60,['p<0.01'],'FontSize', 8)
 elseif P(2)<0.001
     text(0.7,60,['p<0.001'],'FontSize', 8)
 else
     text(0.7,60,['n.s.'],'FontSize', 8)
 end
end
%% 
com2=[par1(idxn) x_weight(idxn)']; 
correlation_matrix(com2,0);


 
%% %%%%%%   LOOK at all cells and pial depth
%% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
    %% EXTACT NECESSARY INFO from structure 
    L23_PC.data= data;
    oris=[0,45,90,135,180,225,270,315];
      field_number = (size(L23_PC.data,1)-1);
        for m=1:field_number
            pial=L23_PC.data{m+1, 3};
            pial_n=pial/max(pial);
            pia_res{:,m}=pial(L23_PC.data{m+1, 4}.resp);
            pia_unres{:,m}=pial(L23_PC.data{m+1, 4}.unres);
            pia_c{:,m}=pial(L23_PC.data{m+1, 4}.contra_only);
            pia_i{:,m}=pial(L23_PC.data{m+1, 4}.ipsi_only);
            pia_b{:,m}=pial(L23_PC.data{m+1, 4}.bino);
            pia_res_n{:,m}=pial_n(L23_PC.data{m+1, 4}.resp);
            pia_unres_n{:,m}=pial_n(L23_PC.data{m+1, 4}.unres);
            pia_c_n{:,m}=pial_n(L23_PC.data{m+1, 4}.contra_only);
            pia_i_n{:,m}=pial_n(L23_PC.data{m+1, 4}.ipsi_only);
            pia_b_n{:,m}=pial_n(L23_PC.data{m+1, 4}.bino);
            pial=[];
            pial_n=[];
            gOSI_c_all{:,m}=L23_PC.data{m+1, 4}.gOSI_c;
            gOSI_i_all{:,m}=L23_PC.data{m+1, 4}.gOSI_i;
            gOSI_b_all{:,m}=L23_PC.data{m+1, 4}.gOSI_b;
            gDSI_c_all{:,m}=L23_PC.data{m+1, 4}.gDSI_c;
            gDSI_i_all{:,m}=L23_PC.data{m+1, 4}.gDSI_i;
            gDSI_b_all{:,m}=L23_PC.data{m+1, 4}.gDSI_b;
            gODI{:,m}=L23_PC.data{m+1, 4}.ODI_r;
            gODI_c_all{:,m}=L23_PC.data{m+1, 4}.ODI_c;
            gODI_i_all{:,m}=L23_PC.data{m+1, 4}.ODI_i;
            gODI_b_all{:,m}=L23_PC.data{m+1, 4}.ODI_b; 
            
            prefOri_ipsi{:,m}=L23_PC.data{m+1, 4}.prefOri_ipsi;
            prefOri_contra{:,m}=L23_PC.data{m+1, 4}.prefOri_contra;
            prefOri_bino_c{:,m}=L23_PC.data{m+1, 4}.prefOri_bino_c;
            prefOri_bino_i{:,m}=L23_PC.data{m+1, 4}.prefOri_bino_i;
            
            prefDir_ipsi{:,m}=L23_PC.data{m+1, 4}.prefDir_ipsi;
            prefDir_contra{:,m}=L23_PC.data{m+1, 4}.prefDir_contra;
            prefDir_bino_c{:,m}=L23_PC.data{m+1, 4}.prefDir_bino_c;
            prefDir_bino_i{:,m}=L23_PC.data{m+1, 4}.prefDir_bino_i;
            
             oppResp_ipsi{:,m}=L23_PC.data{m+1, 4}.oppResp_ipsi;
            oppResp_contra{:,m}=L23_PC.data{m+1, 4}.oppResp_contra;
            oppResp_bino_c{:,m}=L23_PC.data{m+1, 4}.oppResp_bino_c;
            oppResp_bino_i{:,m}=L23_PC.data{m+1, 4}.oppResp_bino_i;
            
            sigma_ipsi{:,m}=L23_PC.data{m+1, 4}.sigma_ipsi;
            sigma_contra{:,m}=L23_PC.data{m+1, 4}.sigma_contra;
            sigma_bino_c{:,m}=L23_PC.data{m+1, 4}.sigma_bino_c;
            sigma_bino_i{:,m}=L23_PC.data{m+1, 4}.sigma_bino_i;
            
            error_ipsi{:,m}=L23_PC.data{m+1, 4}.error_ipsi;
            error_contra{:,m}=L23_PC.data{m+1, 4}.error_contra;
            error_bino_c{:,m}=L23_PC.data{m+1, 4}.error_bino_c;
            error_bino_i{:,m}=L23_PC.data{m+1, 4}.error_bino_i;
            
            rsquare_ipsi{:,m}=L23_PC.data{m+1, 4}.rsquare_ipsi;
            rsquare_contra{:,m}=L23_PC.data{m+1, 4}.rsquare_contra;
            rsquare_bino_c{:,m}=L23_PC.data{m+1, 4}.rsquare_bino_c;
            rsquare_bino_i{:,m}=L23_PC.data{m+1, 4}.rsquare_bino_i;
            
            vari_contra{:,m}=L23_PC.data{m+1, 4}.vari_contra;
            vari_ipsi{:,m}=L23_PC.data{m+1, 4}.vari_ipsi;
            vari_bino_c{:,m}=L23_PC.data{m+1, 4}.vari_bino_c;
            vari_bino_i{:,m}=L23_PC.data{m+1, 4}.vari_bino_i;
            
            varipref_contra{:,m}=L23_PC.data{m+1, 4}.varipref_contra;
            varipref_ipsi{:,m}=L23_PC.data{m+1, 4}.varipref_ipsi;
            varipref_bino_c{:,m}=L23_PC.data{m+1, 4}.varipref_bino_c;
            varipref_bino_i{:,m}=L23_PC.data{m+1, 4}.varipref_bino_i;
            
            ave_c{:,m}=max(L23_PC.data{m+1, 4}.ave_contra);
            ave_i{:,m}=max(L23_PC.data{m+1, 4}.ave_ipsi);
            z_ave_c{:,m}=max(L23_PC.data{m+1, 4}.z_ave_contra);
            z_ave_i{:,m}=max(L23_PC.data{m+1, 4}.z_ave_ipsi);
            ave_bc{:,m}=max(L23_PC.data{m+1, 4}.ave_bino_c);
            ave_bi{:,m}=max(L23_PC.data{m+1, 4}.ave_bino_i);
            z_ave_bc{:,m}=max(L23_PC.data{m+1, 4}.z_ave_bino_c);
            z_ave_bi{:,m}=max(L23_PC.data{m+1, 4}.z_ave_bino_i);
            
            STFT_all_res{:,m}=L23_PC.data{m+1, 5}.res;

        end
%% Discretize pial depth and orientations
%call function discretize_pia_ref
%[ori_pref] = discretize_pia_pref(prefOri_contra,prefOri_ipsi,prefOri_bino_c,prefOri_bino_i,prefDir_contra,prefDir_ipsi,prefDir_bino_c,prefDir_bino_i,pia_c,pia_i,pia_b,1,0,[100:90:400],[0:25:200],[0:40:365]);
[ori_pref] = discretize_pia_pref(gODI_b_all,prefOri_contra,prefOri_ipsi,prefOri_bino_c,prefOri_bino_i,prefDir_contra,prefDir_ipsi,prefDir_bino_c,prefDir_bino_i,pia_c,pia_i,pia_b,1,0,[100:90:400],[0:45:200],[0:40:365]);
%% Discretize pial depth and directions
[dir_pref] = discretize_pia_pref(gODI_b_all,prefOri_contra,prefOri_ipsi,prefOri_bino_c,prefOri_bino_i,prefDir_contra,prefDir_ipsi,prefDir_bino_c,prefDir_bino_i,pia_c,pia_i,pia_b,0,1,[100:90:400],[0:25:180],[0:90:360]);   
%% Plot pref oris vs pial depth 
fig2 = figure;
set(fig2, 'Name', 'Preferred oris');
set(fig2, 'Position', [200, 0, 1500, 1000]);
set(gcf,'color','w');
ori_bin=[0:45:135];
ori_bin=ori_bin-90;
subplot(2,2,1);plot(ori_bin,ori_pref.fract1c,'--o');hold on;plot(ori_bin,ori_pref.fract2c,'--o');hold on;plot(ori_bin,ori_pref.fract3c,'--o');set(gca,'box','off');xlabel('Preferred Orientation');ylabel('Fraction of cells');ylim([0 0.6]);
legend('contra upper L2/3','contra middle L2/3','contra lower L2/3');
subplot(2,2,2);plot(ori_bin,ori_pref.fract1i,'--o');hold on;plot(ori_bin,ori_pref.fract2i,'--o');hold on;plot(ori_bin,ori_pref.fract3i,'--o');set(gca,'box','off');xlabel('Preferred Orientation');ylabel('Fraction of cells');ylim([0 0.6]);
legend('ipsi upper L2/3','ipsi middle L2/3','ipsi lower L2/3');
% subplot(2,2,3);plot(ori_bin,ori_pref.fract1bc,'--o');hold on;plot(ori_bin,ori_pref.fract2bc,'--o');hold on;plot(ori_bin,ori_pref.fract3bc,'--o');set(gca,'box','off');xlabel('Preferred Orientation');ylabel('Fraction of cells');ylim([0 0.4]);
% legend('bino c upper L2/3','bino c middle L2/3','bino c lower L2/3');
subplot(2,2,3);plot(ori_bin,ori_pref.fract1bi,'--o');hold on;plot(ori_bin,ori_pref.fract2bi,'--o');hold on;plot(ori_bin,ori_pref.fract3bi,'--o');set(gca,'box','off');xlabel('Preferred Orientation');ylabel('Fraction of cells');ylim([0 0.6]);
legend('bino upper L2/3','bino middle L2/3','bino lower L2/3');           
figure;scatter(abs(horzcat(prefOri_contra{1,:})),horzcat(pia_c{1,:}));hold on;scatter(abs(horzcat(prefOri_ipsi{1,:})),horzcat(pia_i{1,:}));hold on;scatter(abs(horzcat(prefOri_bino_c{1,:})),horzcat(pia_b{1,:}));
hold on;scatter(abs(horzcat(prefOri_bino_i{1,:})),horzcat(pia_b{1,:}))
xlabel('Preferred Orientation')
set(gca, 'YDir','reverse')
ylabel('Pial distance (µm)');
set(gca,'box','off');
set(gcf,'color','w');
legend('Contra only','Ipsi only','Bino');
%% Add ipsi, contra and bino cells together 
up=[ori_pref.fract1c;ori_pref.fract1i; ori_pref.fract1bi];
mid=[ori_pref.fract2c;ori_pref.fract2i; ori_pref.fract2bi];
down=[ori_pref.fract3c;ori_pref.fract3i; ori_pref.fract3bi];
cellnr=ori_pref.pref_ori_c+ori_pref.pref_ori_i+ori_pref.pref_ori_b;
figure;
set(gcf,'color','w');
errorbar(nanmean(up),nanstd(up)/sqrt(length(up)),'--ro');hold on;errorbar(nanmean(mid),nanstd(mid)/sqrt(length(mid)),'--bo');hold on;errorbar(nanmean(down),nanstd(down)/sqrt(length(down)),'--ko');
xlabel('Preferred Orientation')
ylabel('Fraction of cells (%)');
set(gca,'box','off');
set(gcf,'color','w');
legend('upper L2/3','mid L2/3','lower L2/3'); 
text(1.1,0.5,['N=' num2str(length(prefOri_bino_i))]);
text(1.1,0.48,['n=' num2str(cellnr) ' (c=' num2str(ori_pref.pref_ori_c)  ' i=' num2str(ori_pref.pref_ori_i) ' b=' num2str(ori_pref.pref_ori_b) ')']);
xticks([1:1:5]);
%set(gca,'XTickLabel',{'0-25°' '25-50°','50-75°','75-100°','100-125°','125-150°','150-175°','175-200°'});
%set(gca,'XTickLabel',{'-90 - -65°' '-64 - -40°','-39 - -15°','-14 - 10°','11 - 35°','36 - 60°','61 - 85°','86 - 110°'});
set(gca,'XTickLabel',{'-67.5°' '-22.5°','22.5°','67.5°'});
ylim([0 0.5]);
yticks([0:0.1:0.5]);
set(gca,'YTickLabel',{'0' '10','20','30','40','50'});
legend boxoff;
%% Extract per animal
[ori_prefall] = discretize_pia_pref_peranimal(gODI_b_all,prefOri_contra,prefOri_ipsi,prefOri_bino_c,prefOri_bino_i,prefDir_contra,prefDir_ipsi,prefDir_bino_c,prefDir_bino_i,pia_c,pia_i,pia_b,1,0,[100:90:400],[0:45:200],[0:40:365]);
%% Plot per animal 
up=ori_prefall.fract1_peranimal;
mid=ori_prefall.fract2_peranimal;
down=ori_prefall.fract3_peranimal;
figure;
set(gcf,'color','w');
errorbar(nanmean(up,2),nanstd(up,0,2)/sqrt(length(up)),'-go');hold on;errorbar(nanmean(mid,2),nanstd(mid,0,2)/sqrt(length(mid)),'-bo');hold on;errorbar(nanmean(down,2),nanstd(down,0,2)/sqrt(length(down)),'-ko');
xlabel('Orientation (deg)')
ylabel('Fraction of cells (%)');
set(gca,'box','off');
set(gcf,'color','w');
legend(['upper L2/3' ' (N=' num2str(sum(~isnan(up(1,:)))) ')'],['mid L2/3' ' (N=' num2str(sum(~isnan(mid(1,:)))) ')'],['lower L2/3'' (N=' num2str(sum(~isnan(down(1,:)))) ')']); 
xticks([1:1:5]);
%set(gca,'XTickLabel',{'0-25°' '25-50°','50-75°','75-100°','100-125°','125-150°','150-175°','175-200°'});
%set(gca,'XTickLabel',{'-90 - -65°' '-64 - -40°','-39 - -15°','-14 - 10°','11 - 35°','36 - 60°','61 - 85°','86 - 110°'});
set(gca,'XTickLabel',{'-67.5°' '-22.5°','22.5°','67.5°'});
ylim([0 0.5]);
yticks([0:0.1:0.5]);
set(gca,'YTickLabel',{'0' '10','20','30','40','50'});
legend boxoff;
%% Calculate average across animals independet of pial depth
oripref_all=ori_prefall.bin_ori_idx;
ori_bin=[0:45:200];
for k=1:length(oripref_all)
               temp=oripref_all{:,k};
            
           for l=1:length(ori_bin)-1;
           fract1c(l)=sum(temp==l)/length(temp);
      
           end
           
            fract_ori(:,k)=fract1c;
           end
  %% Barplot of orienatation prference across all mice
figure; bar(nanmean(fract_ori,2))
set(gcf,'color','w');
box off;
xticks([1:1:4]);
set(gca,'XTickLabel',{'-67.5°' '-22.5°','22.5°','67.5°'});
hold on;
err=nanstd(fract_ori,0,2)/sqrt(length(fract_ori));
x=1:4
er = errorbar(x,nanmean(fract_ori,2),err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
xlabel('Preferred Orientation')
ylabel('Fraction of cells (%)');
ylim([0 0.4]);
yticks([0:0.1:0.4]);

set(gca,'YTickLabel',{'0' '10','20','30','40','50'});
text(0.1,0.38,['N=' num2str(length(prefOri_bino_i))]);



%% Direction 
%% Plot pref Dir vs pial depth 
fig2 = figure;
set(fig2, 'Name', 'Preferred Dir');
set(fig2, 'Position', [200, 0, 1500, 1000]);
set(gcf,'color','w');
ori_bin=[0:95:360];
subplot(2,2,1);plot(ori_bin,dir_pref.fract1c,'--o');hold on;plot(ori_bin,dir_pref.fract2c,'--o');hold on;plot(ori_bin,dir_pref.fract3c,'--o');set(gca,'box','off');xlabel('Preferred Direction');ylabel('Fraction of cells');ylim([0 0.6]);
legend('contra upper L2/3','contra middle L2/3','contra lower L2/3');
subplot(2,2,2);plot(ori_bin,dir_pref.fract1i,'--o');hold on;plot(ori_bin,dir_pref.fract2i,'--o');hold on;plot(ori_bin,dir_pref.fract3i,'--o');set(gca,'box','off');xlabel('Preferred Direction');ylabel('Fraction of cells');ylim([0 0.6]);
legend('ipsi upper L2/3','ipsi middle L2/3','ipsi lower L2/3');
% subplot(2,2,3);plot(ori_bin,ori_pref.fract1bc,'--o');hold on;plot(ori_bin,ori_pref.fract2bc,'--o');hold on;plot(ori_bin,ori_pref.fract3bc,'--o');set(gca,'box','off');xlabel('Preferred Orientation');ylabel('Fraction of cells');ylim([0 0.4]);
% legend('bino c upper L2/3','bino c middle L2/3','bino c lower L2/3');
subplot(2,2,3);plot(ori_bin,dir_pref.fract1bi,'--o');hold on;plot(ori_bin,dir_pref.fract2bi,'--o');hold on;plot(ori_bin,dir_pref.fract3bi,'--o');set(gca,'box','off');xlabel('Preferred Direction');ylabel('Fraction of cells');ylim([0 0.6]);
legend('bino upper L2/3','bino middle L2/3','bino lower L2/3');           
figure;scatter(abs(horzcat(prefDir_contra{1,:})),horzcat(pia_c{1,:}));hold on;scatter(abs(horzcat(prefDir_ipsi{1,:})),horzcat(pia_i{1,:}));hold on;scatter(abs(horzcat(prefDir_bino_c{1,:})),horzcat(pia_b{1,:}));
hold on;scatter(abs(horzcat(prefDir_bino_i{1,:})),horzcat(pia_b{1,:}))
xlabel('Preferred Direction')
set(gca, 'YDir','reverse')
ylabel('Pial distance (µm)');
set(gca,'box','off');
set(gcf,'color','w');
legend('Contra only','Ipsi only','Bino');


%% Plot histogram overview of all ODI, OSI, DSI, ori pref, dir pref,
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 400, 800, 300]);
set(gcf,'color','w');
subplot(1,4,1)
h1=histogram(horzcat(gODI{1,:}),'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_ODI(find(~isnan(iv_ODI))),'LineWidth',1);h2.FaceColor = 'b';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);ylabel('Fraction of cells','FontSize',10);xlabel('ODI','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.5;
h2.Normalization = 'probability';
h2.BinWidth = 0.5;
legend(['all (n=' num2str(length(horzcat(gODI{1,:}))) ')'], ['iviv (n=' num2str(sum(~isnan(iv_ODI))) ')']); legend boxoff;

subplot(1,4,3)
h1=histogram([horzcat(gOSI_c_all{1,:}) horzcat(gOSI_i_all{1,:}) horzcat(gOSI_b_all{1,:})],'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_OSI(find(~isnan(iv_OSI))),'LineWidth',1);h2.FaceColor = 'b';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('OSI','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.25;
h2.Normalization = 'probability';
h2.BinWidth = 0.25;

subplot(1,4,4)
h1=histogram([horzcat(gDSI_c_all{1,:}) horzcat(gDSI_i_all{1,:}) horzcat(gDSI_b_all{1,:})],'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_DSI(find(~isnan(iv_DSI))),'LineWidth',1);h2.FaceColor = 'b';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('DSI','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.25;
h2.Normalization = 'probability';
h2.BinWidth = 0.25;

%  
%% CDF plots insteadfigure
% fig22=figure;
% set(fig22, 'Name', 'Binocular protocol');set(fig22, 'Position', [200, 800, 800, 300]);
% set(gcf,'color','w');subplot(1,4,1);cd1=cdfplot(horzcat(gODI{1,:}));cd1.Color='c';hold on; cd2=cdfplot(iv_ODI(find(~isnan(iv_ODI))));cd2.Color='k';
% grid off;box off;title('');ylabel('Cumulative fraction');xlabel('ODI');
% subplot(1,4,2);cd1=cdfplot([horzcat(gOSI_c_all{1,:}) horzcat(gOSI_i_all{1,:}) horzcat(gOSI_b_all{1,:})]);cd1.Color='c';hold on; cd2=cdfplot(iv_OSI(find(~isnan(iv_OSI))));cd2.Color='k';
% grid off;box off;title('');ylabel('Cumulative fraction');xlabel('OSI');

%% Plot histogram overview of all Ori and Dir preference
ODIb=horzcat(gODI_b_all{1,:});
pref_bino_c=horzcat(prefOri_bino_c{1,:});
pref_bino_i=horzcat(prefOri_bino_i{1,:});
pref_bino_c_dir=horzcat(prefDir_bino_c{1,:});
pref_bino_i_dir=horzcat(prefDir_bino_i{1,:});
for m=1:length(ODIb);   
if ODIb>0
    pref_bino(:,m)=pref_bino_c(m);
    pref_bino_dir=pref_bino_c_dir(m);
else ODIb<0
    pref_bino(:,m)=pref_bino_i(m);
     pref_bino_dir=pref_bino_i_dir(m);
end
end
%% Plot histogram overview of all Ori and Dir preference
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 400, 200]);
set(gcf,'color','w');
subplot(1,2,1)
h1=histogram([horzcat(prefOri_ipsi{1,:}) horzcat(prefOri_contra{1,:}) pref_bino]-90,'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(oripref(find(~isnan(oripref)))-90,'LineWidth',1);h2.FaceColor = 'b';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('Orientation (deg)','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 45;
h2.Normalization = 'probability';
h2.BinWidth = 45;
legend(['all (n=' num2str(length(horzcat(gODI{1,:}))) ')'], ['iviv (n=' num2str(sum(~isnan(oripref))) ')']); legend boxoff;
xlim([-90 90]);

subplot(1,2,2);
h1=histogram([horzcat(prefDir_ipsi{1,:}) horzcat(prefDir_contra{1,:}) pref_bino_dir]-180,'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(dirpref(find(~isnan(dirpref)))-180,'LineWidth',1);h2.FaceColor = 'b';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('Direction (deg)','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 90;
h2.Normalization = 'probability';
h2.BinWidth = 90;
xlim([-180 180]);
 
%% ALL Spon activity
%% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% Plot histograms for Popcop and SAD
fig22 = figure;
set(fig22, 'Name', 'Spon activity');set(fig22, 'Position', [200, 200, 400, 200]);
set(gcf,'color','w');

subplot(1,2,1)
h1=histogram(vertcat(popcop_cell{:,1}),'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_popcop(find(~isnan(iv_popcop))),'LineWidth',1);h2.FaceColor = 'm';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('PCI','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.75;
h2.Normalization = 'probability';
h2.BinWidth = 0.75;
legend(['all (n=' num2str(length(vertcat(popcop_cell{:,1}))) ')'], ['iviv (n=' num2str(sum(~isnan(iv_popcop))) ')']); legend boxoff;
ylabel('Fraction of cells','FontSize',10)

subplot(1,2,2)
h1=histogram(horzcat(event_store{:,1}),'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_spon(find(~isnan(iv_spon))),'LineWidth',1);h2.FaceColor = 'm';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('SAD (Hz)','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.03;
h2.Normalization = 'probability';
h2.BinWidth = 0.03;
%% %% ALL SFTF
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file

%% 
SF=[0.02 0.08 0.16];
TF=[1 2 4];
sf_all={};
tf_all={};
for k=1:length(sftf_mat)
    tf_c=[];
    sf_c=[];
    temp=sftf_mat(k).data;
 
for i=1:length(temp)
     if ~isnan(temp(i,1));
       temp1=squeeze(nanmean(max(temp(i,:,:,:),[],4),1));
        [M,I] = max(temp1(:));
     [I_TF, I_SF] = ind2sub(size(temp1),I);
     tf_c(:,i)=TF(I_TF);
     sf_c(:,i)=SF(I_SF);
     else
        tf_c(:,i)=NaN;
        sf_c(:,i)=NaN;
     end
end

temp=[];
sf_all{k}=sf_c;
tf_all{k}=tf_c;
end
%% 

%% 

fig23 = figure;
set(fig23, 'Name', 'SFTF');set(fig23, 'Position', [200, 800, 400, 200]);
set(gcf,'color','w');
subplot(1,2,1);imagesc(sftfall);c=colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);caxis([0 0.25]);
set(gca,'YTickLabel',{'1' '2','4'});set(gca,'Ydir','reverse');ylabel('TF');
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;xtickangle(45);
hold on; title(['all' ' n=' num2str(length(find(~isnan(sf_com))))]);caxis([0 0.25]);
subplot(1,2,2);imagesc(sftfalli);colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);
set(gca,'YTickLabel',{'1' '2','4'});set(gca,'Ydir','reverse');ylabel('TF');
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;
hold on; title(['iviv' ' n=' num2str(length(find(~isnan(s_SF))))]);caxis([0 0.25]);
xtickangle(45);
%% histogram


fig23 = figure;
set(fig23, 'Name', 'SFTF');set(fig23, 'Position', [200, 400, 400, 200]);
set(gcf,'color','w');

subplot(1,2,1)
h1=histogram(sf_com,'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_SF(find(~isnan(iv_SF))),'LineWidth',1);h2.FaceColor = 'm';h2.FaceAlpha =0.25;box off;ylim([0 8]);xlabel('SF','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.75;
h2.Normalization = 'probability';
h2.BinWidth = 0.75;
legend(['all (n=' num2str(length(vertcat(popcop_cell{:,1}))) ')'], ['iviv (n=' num2str(sum(~isnan(iv_popcop))) ')']); legend boxoff;
ylabel('Fraction of cells','FontSize',10)

subplot(1,2,2)
h1=histogram(horzcat(event_store{:,1}),'LineWidth',1);h1.FaceColor = 'w';box off;hold on;h2=histogram(iv_spon(find(~isnan(iv_spon))),'LineWidth',1);h2.FaceColor = 'm';h2.FaceAlpha =0.25;box off;ylim([0 0.6]);xlabel('SAD (Hz)','FontSize',10);
h1.Normalization = 'probability';
h1.BinWidth = 0.03;
h2.Normalization = 'probability';
h2.BinWidth = 0.03;







%% 
%% Morpho density corr basal
% calculate the correlation between morpho density and maps
morphoMaps = {str(non_nan_cells).morphoMap_basal};
%morphoMaps = {str(non_nan_cells).morphoMap_apical};
excMaps = {str(non_nan_cells).excMap};
inhMaps = {str(non_nan_cells).inhMap};
non_nan_names = {str(non_nan_cells).cellName};
correlation_values = cell(cell_num,3);
% for all the cells
for cells = 1:cell_num
    if isempty(morphoMaps{cells})
        correlation_values{cells,3} = 'nanCell';
        continue
    end
    % save the name of the cell
    correlation_values{cells,3} = non_nan_names{cells};
    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = excMaps{cells}*-1;
            case 2
                map = inhMaps{cells};
        end
        % calculate the correlation and store
        correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
        %correlation_values_ap{cells,maps} = corr(map(:), morphoMaps{cells}(:));
    end  
end

%% Morpho density corr apical
% calculate the correlation between morpho density and maps
morphoMaps = {str(non_nan_cells).morphoMap_apical};
excMaps = {str(non_nan_cells).excMap};
inhMaps = {str(non_nan_cells).inhMap};
non_nan_names = {str(non_nan_cells).cellName};
correlation_values_ap = cell(cell_num,3);
% for all the cells
for cells = 1:cell_num
    if isempty(morphoMaps{cells})
        correlation_values{cells,3} = 'nanCell';
        continue
    end
    % save the name of the cell
    correlation_values{cells,3} = non_nan_names{cells};
    % for the 2 types of map
    for maps = 1:2
        % get the map
        switch maps
            case 1
                map = excMaps{cells}*-1;
            case 2
                map = inhMaps{cells};
        end
        % calculate the correlation and store
        %correlation_values{cells,maps} = corr(map(:), morphoMaps{cells}(:));
        correlation_values_ap{cells,maps} = corr(map(:), morphoMaps{cells}(:));
    end  
end
%% 

figure;
a=81
       m=plot_tree(str(a).morphtraces{1,1},[1 0 0],[0 str(a).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(a).morphtraces{1,2},[0 0 0],[0 str(a).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(a).morphtraces{1,3},[0 0 1],[0 str(a).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
       
set(gca,'Ydir','reverse')
%% 

test=find(iviv_cells==1)'
test2=find(morpho_cells==1)'
[mo mo2]=intersect(test,test2)

fig23=figure;
set(fig23, 'Name', 'Spon activity');set(fig23, 'Position', [200, 200, 800, 800]);hold on;set(gcf,'color','w');
for a=1:length(mo)
subplot(6,6,a)
 m=plot_tree(str(mo(a)).morphtraces{1,1},[1 0 0],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(mo(a)).morphtraces{1,2},[0 0 0],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(mo(a)).morphtraces{1,3},[0 0 1],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
set(gca,'Ydir','reverse');
end
%% 
fig23=figure;
set(fig23, 'Name', 'Spon activity');set(fig23, 'Position', [200, 200, 800, 800]);hold on;set(gcf,'color','w');
for a=1:length(mo)
subplot(6,6,a)
 imagesc(str(mo(a)).morphoMap_apical);%hold on;imagesc(str(mo(a)).morphoMap_basal);
set(gca,'Ydir','reverse');
title(num2str(str(mo(a)).sliceOri))
hh(a)=str(mo(a)).sliceOri;
end

fig23=figure;
set(fig23, 'Name', 'Spon activity');set(fig23, 'Position', [200, 200, 800, 800]);hold on;set(gcf,'color','w');
for a=1:length(mo)
subplot(6,6,a)
 imagesc(str(mo(a)).morphoMap_basal);%hold on;imagesc(str(mo(a)).morphoMap_basal);
set(gca,'Ydir','reverse');
title(num2str(str(mo(a)).cellName))
end

%% 
 [a1 a2] = intersect(non_nan_idx,test2);
 [b1 b2] = intersect(non_nan_idx,test);
 
 [c1 c2] = intersect(a1,b1);
 
 [d1 d2]= intersect(non_nan_idx,a1)
 
 %% Correlation morphology and input/iviv
 

for i=1:length(a1)
    
        LMO(i) = sum(sum(str(a1(i)).morphoMap_basal(:,4:8)));
        RMO(i) = sum(sum(str(a1(i)).morphoMap_basal(:,9:13)));
        LMO_a(i) = sum(sum(str(a1(i)).morphoMap_apical(:,4:8)));
        RMO_a(i) = sum(sum(str(a1(i)).morphoMap_apical(:,9:13)));      
       vr(i) = str(a1(i)).sliceOri;   
end
rm_lm=(RMO-LMO)./(RMO+LMO)
rm_lm_a=(RMO_a-LMO_a)./(RMO_a+LMO_a);

%% 

for i=1:length(vr)
    if vr(i)==0
    rm_lm_f(i)=rm_lm(i)*-1;
    rm_lm_fa(i)=rm_lm_a(i)*-1
    else
    rm_lm_f(i)=rm_lm(i)*1  
    rm_lm_fa(i)=rm_lm_a(i)*1  
    end
end
%% 
pc2in=cell_cell(d2,5);
figure;scatter(rm_lm_f,pc2in);

[R P]=corrcoef(rm_lm_f,pc2in,'row','complete');

%% 

 %% Correlation morphology and input/iviv
 
LMO = zeros(length(non_nan_idx),1);
RMO = zeros(length(non_nan_idx),1);
for i=1:length(non_nan_idx)
    
    if ~isempty(str(non_nan_idx(i)).morph)==1
%         LMO(i)=sum(sum(str(non_nan_idx(i)).morphoMap_basal(4:6,4:7)));
%         RMO(i)=sum(sum(str(non_nan_idx(i)).morphoMap_basal(4:6,9:12)));
        LMO(i) = sum(sum(str(non_nan_idx(i)).morphoMap_basal(:,4:8)));
        RMO(i) = sum(sum(str(non_nan_idx(i)).morphoMap_basal(:,9:13)));
        LMO_a(i) = sum(sum(str(non_nan_idx(i)).morphoMap_apical(:,4:8)));
        RMO_a(i) = sum(sum(str(non_nan_idx(i)).morphoMap_apical(:,9:13)));
     
    else
        LMO(i)=NaN;
        RMO(i)=NaN;
        LMO_a(i)=NaN;
        RMO_a(i)=NaN;
    end
end
rm_lm=(RMO-LMO)./(RMO+LMO)
rm_lm_a=(RMO_a-LMO_a)./(RMO_a+LMO_a);
%% Offset of apical and basal dendrite 
%BASAL DENDRITE 
fig21= figure;
set(fig21, 'Name', 'Offset apical and basal');set(fig21, 'Position', [200, 0, 1000, 400]);set(gcf,'color','w');
par1=rm_lm;
par2=oripref-90;
%par2=iv_DSI;
%par3=pia_input;
par3=nan_vector;
[cmap]=buildcmap('mwg');
pointsize=20;
subplot(1,2,1);
scatter(par1,par2,pointsize,par3,'filled');
set(gcf,'color','w');
box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth'; 
idxn = find(~isnan(par2));
temp1=par1(idxn);
temp2=par2(idxn)';
idxn2 = find(~isnan(temp1));
P = polyfit(temp1(idxn2),temp2(idxn2),1);
    yfit = P(1)*temp1(idxn2)+P(2);
    hold on;box off;
    plot(temp1(idxn2),yfit,'--k');set(gca,'box','off');set(gcf,'color','w');
 ylabel('Ori pref(°)');xlabel('Offset basal dendrite');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(0.75,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(0.75,70,['p<0.05'])

 elseif P(2)<0.01 & P(2)>0.001
      text(0.75,70,['p<0.01'])
 elseif P(2)<0.001
     text(0.75,70,['p<0.001'])
 else
     text(0.75,70,['n.s.'])
 end

  %% to check cell IDS
 figure;scatter3(par1,par2,par3,'filled');
 %% 
 
%APICAL DENDRITE 
subplot(1,2,2);
par1=[];
par2=[];
par2=oripref-90;
par1=rm_lm_a;
[cmap]=buildcmap('mwg');
pointsize=20;
scatter(par1,par2,pointsize,par3,'filled');
set(gcf,'color','w');
box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth'; 
idxn = find(~isnan(par2));
temp1=par1(idxn);
temp2=par2(idxn)';
idxn2 = find(~isnan(temp1));
P = polyfit(temp1(idxn2),temp2(idxn2)',1);
    yfit = P(1)*temp1(idxn2)+P(2);
    hold on;box off;
    plot(temp1(idxn2),yfit,'k');set(gca,'box','off');set(gcf,'color','w');
 ylabel('Ori pref(°)');xlabel('Offset apical dendrite');
 [R P]=corrcoef(par1,par2,'row','complete');
 text(0.75,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(0.75,70,['p<0.05'])

 elseif P(2)<0.01 & P(2)>0.001
      text(0.75,70,['p<0.01'])
 elseif P(2)<0.001
     text(0.75,70,['p<0.001'])
 else
     text(0.75,70,['n.s.'])
 end
 %% 
  %% Morpho input correlation to iv

 for i=1:length(morphoMaps)
 if ~isempty(morphoMaps{i})==1;
     morpht(i)=1;
     exmi(i)=correlation_values{i,1};
     inmi(i)=correlation_values{i,2};
     exmi_ap(i)=correlation_values_ap{i,1};
     inmi_ap(i)=correlation_values_ap{i,2};
 else isempty(morphoMaps{i})==1;
     morpht(i)=NaN;
     exmi(i)=NaN;
     inmi(i)=NaN;
     exmi_ap(i)=NaN;
     inmi_ap(i)=NaN;
 end
 end
 %% 
par1=cell_cell(:,1);
%par1=
%par2=oripref-90;
par2=exmi;
par3=oripref-90;
[cmap]=buildcmap('mwg');
pointsize=20;
figure;scatter(par1,par2,pointsize,par3,'filled');
 set(gcf,'color','w');
 box off; c=colorbar;colormap(cmap);c.Label.String = 'Pial depth'; 
idxn = find(~isnan(par2));
P = polyfit(par1(idxn),par2(idxn)',1);
    yfit = P(1)*par1+P(2);
    hold on;box off;
    plot(par1,yfit,'k');set(gca,'box','off');set(gcf,'color','w');axis square;
 ylabel('Ori pref(°)');xlabel('PC2in');
 
 [R P]=corrcoef(par1,par2,'row','complete');
 text(1.5,80,['r= ' mat2str(round(R(2),2))]);
 if P(2)<0.05 & P(2)>0.01
     text(1.5,70,['p<0.05'])

 elseif P(2)<0.01 & P(2)>0.001
      text(1.5,70,['p<0.01'])
 elseif P(2)<0.001
     text(1.5,70,['p<0.001'])
 else
     text(1.5,70,['n.s.'])
 end

 
 %% 
 mo=m_id;
 fig23=figure;
set(fig23, 'Name', 'Spon activity');set(fig23, 'Position', [200, 200, 800, 800]);hold on;set(gcf,'color','w');
for a=1:length(mo)
subplot(8,8,a)
 m=plot_tree(str(mo(a)).morphtraces{1,1},[1 0 0],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'r'
       m1=plot_tree(str(mo(a)).morphtraces{1,2},[0 0 0],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'k'
       m2=plot_tree(str(mo(a)).morphtraces{1,3},[0 0 1],[0 str(mo(a)).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'b'
      hold on;
      title({str(mo(a)).cellName}')
set(gca,'Ydir','reverse');
end
%% 
fig23=figure;
set(fig23, 'Name', 'Spon activity');set(fig23, 'Position', [200, 200, 800, 800]);hold on;set(gcf,'color','w');
for a=1:length(mo)
subplot(8,8,a)
 imagesc(str(mo(a)).morphoMap_apical);%hold on;imagesc(str(mo(a)).morphoMap_basal);
set(gca,'Ydir','reverse');
title(num2str(str(mo(a)).sliceOri))
hh(a)=str(mo(a)).sliceOri;
end
%% 
for i=1:157
str_m(i).morph_flip_again=m_flip_a(i)
end
%% Load raw map and plot
ccr=find(nan_vector==105);
flipo=slice_ori(ccr);

pathName='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Cellsfor_iviv_paper\Cellsfor_iviv_paper\170613\SW0002\map02'
list=dir([char(pathName) filesep '*.xsg']);
j=1;
load([char(pathName) filesep list(j).name],'-mat');
sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
samples_per_sweep = header.ephys.ephys.traceLength*sr;
timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
traces=data.ephys.trace_1;%raw ephys trace
%Filter
cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter';
traces = lowpassfilt(traces, order, cutoff, sr, type);
ind_traces=reshape(traces,[length(traces)/256 256]);
%Baseline subtraction
base_start=1;
base_end=100;
baseline=ind_traces(base_start*srF:base_end*srF,:);
bs_traces=ind_traces-mean(baseline);%subtract baseline
flipFlag = 0;
flipx = 1;
array=bs_traces;
traceLength = size(array,1 );
mapPattern=header.mapper.mapper.mapPatternArray;
if flipx ==1
    mapPattern=fliplr(mapPattern)
end
 for n=1:numel(mapPattern)
 newArray(:,find(mapPattern==n)) = array(:,n);
 end
startTime = 1;
stopTime = traceLength;
showStart = 99*srF;
showStop = 700*srF;
array = newArray(startTime:stopTime, :);
[rows,cols] = size(array);
totalTime = (rows-1)/sr; 
xTimeAxis = linspace(0, totalTime, rows)';
traceAxis = ( 1 : cols );
[sizeX, sizeY] = size(mapPattern);
yFactor = 500; % offset, in pA
%xFactor = 500
scaleBarTxt = 'pA';
if flipFlag == 1
yFactor = -yFactor; 
end
offsetVector = yFactor * ( 0 : cols-1 );
offsetArray = repmat(offsetVector, rows, 1);
array = array-offsetArray;      
% set up the figure -------------------------------------------------------------
x = .14; 
y = .14; 
w = .8; 
h = .8; 
F=figure;
subplotRows = 1; subplotCols = sizeY; plotnum = 0;
for N = 1:sizeY
startInd = (N-1)*sizeX + 1;
endInd = N*sizeX;
plotnum = plotnum+1;

    %     hsub(plotnum) = subplot(subplotRows,subplotCols,plotnum);
    
%     pos1 = 0.025 + (N - 1)*(0.96/sizeY);
%      pos2 = 0.02;
%     pos3 = 0.05;
%     pos4 = 0.96;

   pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    pos2 = 0.02;
    pos3 = 0.038;
    pos4 = 0.96

    hsub(N) = axes('Position', [pos1 pos2 pos3 pos4]);

%     set(gca, 'Position', );

    plot(xTimeAxis(showStart:showStop), array(showStart:showStop,startInd:endInd),'Color','k');
    hold on;
%                 y1=get(gca,'ylim');
%                 x1= redpeak_start/(srF*100);
%                 hold on;
%                 p1=plot([x1 x1],y1,'--','Color','r');
%                 p1.Color(4) = 0.3;
%                 hold on;
% %                 y1=get(gca,'ylim');
% %                 x1=redpeak_end/(srF*100);
% %                 hold on;
% %                 p2=plot([x1 x1],y1,'--','Color','r');
% %                 p2.Color(4) = 0.3;
% %                 hold on;
%                 %blue vertical lines
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_start/(srF*100);
%                 hold on;
%                 p3=plot([x1 x1],y1,'--','Color','b');
%                 p3.Color(4) = 0.3;
%                 hold on;
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_end/(srF*100);
%                 hold on;
%                 p4=plot([x1 x1],y1,'--','Color','b');
%                 p4.Color(4) = 0.3;
    
minval = min(mean(array(1:100,startInd:endInd)));
maxval = max(mean(array(1:100,startInd:endInd)));
tweakFactor = abs(maxval - minval)*0.05;
yrange = [minval-tweakFactor maxval+tweakFactor];
set(gca, 'YLim', yrange);
set(gca, 'XLim', [(showStart-200)/sr (showStop+200)/sr]);
xlabel('Seconds');
set(gcf,'color','w');
%  if flipo==0
%      set(gca, 'XDir','reverse')
% end   
                
end
          
set(findobj(gcf, 'Type', 'axes'), 'Visible', 'off');
    
% scalebar lines
    Y = mean(array(:,end))+yFactor/4;
% Y = min(get(gca, 'YLim'));
    hscalebar = line([.1 .2], [Y Y]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');
    hscalebar = line([.1 .1], [Y Y+yFactor/2]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');

% scalebar text
    ht(1) = text(.12, Y+yFactor/6, '100 ms'); 
    ht(2) = text(.12, Y+yFactor/3, [num2str(yFactor/2) ' ' scaleBarTxt]); 
    set(ht, 'Color', 'k', 'FontSize', 8, 'Tag', 'scaleBarText');
%axis square;
%% Plot exmple cells in vivo values for figure 1
oris=[0,45,90,135,180,225,270,315];
u=137;
F=figure;
set(gcf, 'Position', [300,600, 800, 200]);
set(gcf,'color','w');
H=subplot(1,4,1);
set(H, 'Position', [0.1, 0.4, 0.1, 0.3]);
if isempty(str(u).morphtraces)==0
m=plot_tree(str(u).morphtraces{1,1},[1 0 0],[0 str(u).morphtraces{1,5} 0],[],1,'-b');
m.EdgeColor = 'r'
hold on;
m1=plot_tree(str(u).morphtraces{1,2},[0 0 0],[0 str(u).morphtraces{1,5} 0],[],1,'-b');hold on;
m1.EdgeColor = 'k'
hold on;
m2=plot_tree(str(u).morphtraces{1,3},[0 0 1],[0 str(u).morphtraces{1,5} 0],[],1,'-b');hold on;
m2.EdgeColor = 'b'
set(gca,'Ydir','reverse');
xlim([-250 250]);
ylim([0 350]);
%axis off;
else
 m=plot(str(u).pialD,'Marker','^','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',4);
 axis off;
end

    
%-----------Polarplot_ipsi_contra
hold on;
H2=subplot(1,4,2);
set(H2, 'Position', [0.3, 0.4, 0.1, 0.3]);
polarplot(deg2rad(oris([1:end 1])),str(u).OD_PSTH_norm_c([1:end 1]));
hold on;
polarplot(deg2rad(oris([1:end 1])),str(u).OD_PSTH_norm_i([1:end 1]));
ax = gca; % current axes
%ax.ThetaTickLabel = [];
rticks([0 0.5 1]);
ax.RColor = [0 0 0];
ax.LineWidth = 2;

    
%-----------SFTF
temp=str(u).SFTF_integ_dir;
temp1=squeeze(nanmean(max(temp(:,:,:,:),[],4),1));
temp2=temp1./max(max(temp1))
hold on;
H3=subplot(1,4,3);
set(H3, 'Position', [0.5, 0.4, 0.1, 0.3]);
imagesc(temp2');c=colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);caxis([0 1]);
 set(gca,'YTickLabel',{'1' '2','4'});%set(gca,'Ydir','reverse');
 ylabel('TF');
 set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;xtickangle(45);

xtickangle(45);


%-----------------spon activity-------------
hold on
H4=subplot(1,4,4);
set(H4, 'Position', [0.7, 0.4, 0.05, 0.3]);
colors = {'k','g'};
x=1:7; [ax,h1,h2]=plotyy(x,str(u).iv_spon,x,log(str(u).iv_popcop), @(x,y) bar(2,str(u).iv_spon,colors{1}), @(x,y) bar(4,log(str(u).iv_popcop),colors{2}));
hold on;
set(ax(1),'YLim',[0 0.2]);
% set(ax(2),'YTicks',[0 300]);
set(ax(2),'YLim',[-2 1]);
% set(ax(1),'YTicks',[0 3]);
% ylabel(ax(2), 'PCI');
% ylabel(ax(1), 'SEC');
set(gca,'Xtick',2:2:4);
set(gca,'XTickLabel',{'SAD','PCI'});
set(ax(1),'ytick',[0 0.2]);
set(ax(2),'ytick',[-2 1]);
box off;
xtickangle(45);
%% Spontaneuous activity plot
%Load single exp adn then do this: exp14395
temp=thres_map(:,3);temp(find(temp==0))=NaN;

figure;set(gcf,'color','w');plot(data(2500:3000,3),'k');hold on;plot(temp(2500:3000)*21,'ro');box off;%axis off;
axis off;
%scale barx
hold on;x1=450;x2=500;p1=plot([x1 x2],[-5 -5],'-','Color','k','LineWidth',2);
%scale bary
hold on;y2=-3 ;y1=-5;p1=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',2); 
%% 
aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
    if ivivc.icontra(i)==1 & ivivc.iipsi(i)==0
        str_m(aiviv(i)).contra_o=1;
    elseif ivivc.icontra(i)==0 & ivivc.iipsi(i)==1
   str_m(aiviv(i)).ipsi_o=1;
    elseif ivivc.icontra(i)==1 & ivivc.iipsi(i)==1
     str_m(aiviv(i)).bino_o=1;   
    else ivivc.icontra(i)==0 & ivivc.iipsi(i)==0
        str_m(aiviv(i)).unres_o=1;
       
    end
end

%% 

aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
   str_m(aiviv(i)).ori_a=[ivivc.iOric(i) ivivc.iOrii(i)];
end
%% 

for i=1:length(aiviv)
   str_m(aiviv(i)).dir_a=[ivivc.iDirc(i) ivivc.iDiri(i)];
end







%% 
%% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 for i=1:length(L23_PC)
     con_all{i,:}=L23_PC(i).OD.contra;
     ips_all{i,:}=L23_PC(i).OD.ipsi;
     ODI_all{i,:}=L23_PC(i).OD.ODI;
     prefOri_all{i,:}=L23_PC(i).OD.prefOri(:,:);
     prefDir_all{i,:}=L23_PC(i).OD.prefDir(:,:);
     OSIall{i,:}=L23_PC(i).OD.gOSI(:,:); 
     DSIall{i,:}=L23_PC(i).OD.gDSI(:,:);   
     sf_all{i,:}=L23_PC(i).SFTF.sf;
     tf_all{i,:}=L23_PC(i).SFTF.tf;
     sad_all{i,:}=L23_PC(i).spon.sad;
     pci_all{i,:}=L23_PC(i).spon.pci;
     sftf_resp_all{i,:}=L23_PC(i).SFTF.ov_resp;
     if isempty(L23_PC(i).SFTF.Fit)==0
     sftf_ori{i,:}=[L23_PC(i).SFTF.Fit.PrefRsp];
     else
     sftf_ori{i,:}=NaN;
     end
 end
 %% 
  con=horzcat(con_all{:});
  ips=horzcat(ips_all{:});
  odi=horzcat(ODI_all{:});
  pOSI_all=vertcat(prefOri_all{:});
  pDSI_all=vertcat(prefDir_all{:});
  OSI_all=vertcat(OSIall{:});
  DSI_all=vertcat(OSIall{:});
  sf_a=horzcat(sf_all{:});
  tf_a=horzcat(tf_all{:});
  sad_a=horzcat(sad_all{:});
  pci_a=vertcat(pci_all{:});
  sftfs_res_a=horzcat(sftf_resp_all{:});
  sftf_ori_a=horzcat(sftf_ori{:});
  %% OD protocol
 
  for i=1:length(odi)
      if con(i)==1 & ips(i)==0;
          
          pORI_a(i)=pOSI_all(i,1);
          pDIR_a(i)=pDSI_all(i,1);
          pOSI_a(i)=OSI_all(i,1);
          pDSI_a(i)=DSI_all(i,1);
          
      elseif con(i)==0 & ips(i)==1;
          pORI_a(i)=pOSI_all(i,2);
          pDIR_a(i)=pDSI_all(i,2);
          pOSI_a(i)=OSI_all(i,2);
          pDSI_a(i)=DSI_all(i,2);
      elseif con(i)==1 & ips(i)==1;
          if odi(i)>0
              pORI_a(i)=pOSI_all(i,1);
              pDIR_a(i)=pDSI_all(i,1);
              pOSI_a(i)=OSI_all(i,1);
              pDSI_a(i)=DSI_all(i,1);
          else odi(i)<0
              pORI_a(i)=pOSI_all(i,2);
              pDIR_a(i)=pDSI_all(i,2);
              pDSI_a(i)=DSI_all(i,2);
              pOSI_a(i)=OSI_all(i,2);
          end
      else con(i)==0 & ips(i)==0
       pORI_a(i)=NaN;
       pDIR_a(i)=NaN;
       pOSI_a(i)=NaN;
       pDSI_a(i)=NaN;
      end
  end
    %% SFTF protocol
 
  for i=1:length(sftfs_res_a)
      if sftfs_res_a(i)==1;
          sf_com(i)=sf_a(i);
          tf_com(i)=tf_a(i);    
      else 
       sf_com(i)=NaN;
       tf_com(i)=NaN;
      end
  end
  


%% spon activity+orientation preference
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 400, 100]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(pOSI_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(iv_OSI(find(~isnan(iv_OSI))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('OSI','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
subplot(1,2,2);
binRange = -90:22.5:90;
hcx = histcounts(pORI_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(oripref(find(~isnan(oripref))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Orientation (deg)','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
%legend(['all (n=' num2str(length(pORI_a(find(~isnan(pORI_a)))))], ['iviv (n=' num2str(sum(~isnan(oripref))) ')']); legend boxoff;
%% 
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 400, 100]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(pDSI_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(iv_DSI(find(~isnan(iv_DSI))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('DSI','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
subplot(1,2,2);
binRange = 0:45:360;
hcx = histcounts(pDIR_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(dirpref(find(~isnan(dirpref))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Direction (deg)','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
xticks([0:45:360]);
xtickangle(45)
%legend(['all (n=' num2str(length(pORI_a(find(~isnan(pORI_a)))))], ['iviv (n=' num2str(sum(~isnan(oripref))) ')']); legend boxof
%% 
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 400, 100]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0.0013:0.025:0.2;
hcx = histcounts(sad_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(iv_spon(find(~isnan(iv_spon))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('SAD (Hz)','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.8]);yticks([0:0.4:0.8]);
subplot(1,2,2);
binRange = -5:0.7:5;
hcx = histcounts(pci_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(iv_popcop(find(~isnan(iv_popcop))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('PCI','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.8]);yticks([0:0.4:0.8]);
%% 
figure;set(gcf,'color','w');
resp_all=contra_cells+ipsi_cells+bino_cells;
b=bar([1 2],[contra_cells/resp_all ipsi_cells/resp_all bino_cells/resp_all;  nansum(iv_contra)/sum(~isnan(iv_ODI)) ...
    nansum(iv_ipsi)/sum(~isnan(iv_ODI)) nansum(iv_bino)/sum(~isnan(iv_ODI))],'stacked');
box off;
 b(1, 1).FaceColor='b';b(1, 2).FaceColor='r';b(1, 3).FaceColor='w';
 xticklabels({'all','iviv'});
 legend('contra','ipsi','bino');legend boxoff
%% SF TF

stf1=length(find(sf_com(find(tf_com==1))==0.16))/length(find(~isnan(sf_com)));
stf2=length(find(sf_com(find(tf_com==2))==0.16))/length(find(~isnan(sf_com)));
stf3=length(find(sf_com(find(tf_com==4))==0.16))/length(find(~isnan(sf_com)));

stf4=length(find(sf_com(find(tf_com==1))==0.08))/length(find(~isnan(sf_com)));
stf5=length(find(sf_com(find(tf_com==2))==0.08))/length(find(~isnan(sf_com)));
stf6=length(find(sf_com(find(tf_com==4))==0.08))/length(find(~isnan(sf_com)));

stf7=length(find(sf_com(find(tf_com==1))==0.02))/length(find(~isnan(sf_com)));
stf8=length(find(sf_com(find(tf_com==2))==0.02))/length(find(~isnan(sf_com)));
stf9=length(find(sf_com(find(tf_com==4))==0.02))/length(find(~isnan(sf_com)));
sftfall=[stf7 stf4 stf1;stf8 stf5 stf2;stf9 stf6 stf3];

s_SF=iv_SF((~isnan(iv_SF)));
s_TF=iv_TF((~isnan(iv_TF)));
stf1i=length(find(s_SF(find(s_TF==1))==0.16))/length(find(~isnan(s_SF)));
stf2i=length(find(s_SF(find(s_TF==2))==0.16))/length(find(~isnan(s_SF)));
stf3i=length(find(s_SF(find(s_TF==4))==0.16))/length(find(~isnan(s_SF)));

stf4i=length(find(s_SF(find(s_TF==1))==0.08))/length(find(~isnan(s_SF)));
stf5i=length(find(s_SF(find(s_TF==2))==0.08))/length(find(~isnan(s_SF)));
stf6i=length(find(s_SF(find(s_TF==4))==0.08))/length(find(~isnan(s_SF)));

stf7i=length(find(s_SF(find(s_TF==1))==0.02))/length(find(~isnan(s_SF)));
stf8i=length(find(s_SF(find(s_TF==2))==0.02))/length(find(~isnan(s_SF)));
stf9i=length(find(s_SF(find(s_TF==4))==0.02))/length(find(~isnan(s_SF)));
sftfalli=[stf7i stf4i stf1i;stf8i stf5i stf2i;stf9i stf6i stf3i];
%% 

fig23 = figure;
set(fig23, 'Name', 'SFTF');set(fig23, 'Position', [200, 800, 400, 200]);
set(gcf,'color','w');
subplot(1,2,1);imagesc(sftfall);c=colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);caxis([0 0.3]);
set(gca,'YTickLabel',{'1' '2','4'});ylabel('TF');%set(gca,'Ydir','reverse')
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;xtickangle(45);
hold on; title(['all' ' n=' num2str(length(find(~isnan(sf_com))))]);caxis([0 0.3]);

subplot(1,2,2);imagesc(sftfalli);colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);
set(gca,'YTickLabel',{'1' '2','4'});ylabel('TF');%set(gca,'Ydir','reverse')
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;
hold on; title(['iviv' ' n=' num2str(length(find(~isnan(s_SF))))]);caxis([0 0.3]);
xtickangle(45);
%% 


fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 400, 100]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0.0013:0.025:0.2;
hcx = histcounts(sad_a,[binRange Inf],'Normalization','probability');
hcy = histcounts(iv_spon(find(~isnan(iv_spon))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('SAD (Hz)','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.8]);yticks([0:0.4:0.8]);
subplot(1,2,2);
binRange = -90:22.5:90;
hcx = histcounts(pOSI_all-90,[binRange Inf],'Normalization','probability');
hcy = histcounts(oripref(find(~isnan(oripref)))-90,[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Orientation (deg)','FontSize',7);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.3]);yticks([0:0.15:0.3]);

%% 

figure;set(gcf,'color','w')

h1=histogram(oripref-90);h1.FaceColor = 'k';box off;hold on;h2=histogram(oripref_sftf-90,'LineWidth',1);h2.FaceColor = 'm';h2.FaceAlpha =0.25;box off;ylim([0 30]);xlabel('Ori pref','FontSize',10);

%legend(['all (n=' num2str(length(vertcat(popcop_cell{:,1}))) ')'], ['iviv (n=' num2str(sum(~isnan(iv_popcop))) ')']); legend boxoff;
ylabel('Cells','FontSize',10)
legend('mono stimulation','bino stimulation')








%% 




figure;hold on;
for i=1:28
temp=ty(:,:,m);
temp2=fracinh(m,:);
temp3=orist(m);
%imagesc(islocalmax(temp(:,:,i),1));
plot(temp2(i,:),'k');
% hold on;
% title(num2str(temp3(i)))
end
hold on
for i=29:54
temp=ty(:,:,m);
temp2=fracinh(m,:);
temp3=orist(m);
%imagesc(islocalmax(temp(:,:,i),1));
plot(temp2(i,:),'m');
% hold on;
% title(num2str(temp3(i)))
end
%% 
figure;hold on;
for i=1:length(idxm1)
    plot(fracinh(idxm1(i),:),'k')
end
hold on;
for i=1:length(idxm1)
    plot(fracinh(idxm4(i),:),'m')
end
