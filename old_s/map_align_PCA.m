function [ aligned_maps_ex aligned_maps_in coeff_ex,score_ex,coeff_in,score_in,coeff_com,score_com] = map_align_PCA(str,c_id)
%% Define array names and get the arrays without NaNs etc.
field_list = {'subpixel_raw_excMap','subpixel_raw_inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
%iviv_cells = [str(:).iviv]==1;
% iviv_cells = [str(:).resp]==1;
%  dsi_cells=[str(:).DSIpref]>0.25;
% morpho_cells = ~cellfun(@isempty, {str.morph});
% cell_idx = find(iviv_cells&morpho_cells);
 cell_idx = c_id;
 %cell_idx=find(dsi_cells);
%cell_idx = 1:length(str);
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
        if contains(field_list{fields}, 'subpixel_raw_excMap')
            field_cell{fields} = field_cell{fields}./min(field_cell{fields});
        elseif contains(field_list{fields}, 'subpixel_raw_inhMap')
            field_cell{fields} = field_cell{fields}./max(field_cell{fields});
        end       
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});    
    
    for fields = 1:field_number
        field_cell_raw{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'subpixel_raw_excMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        elseif contains(field_list{fields}, 'subpixel_raw_inhMap')
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
%% %% Get 16x16 maps for ex and in 
ex_map = reshape(original_maps(:,1:256)',16,16,length(nan_vector));
in_map = reshape(original_maps(:,257:512)',16,16,length(nan_vector));
% Get 16x16 maps for ex and in RAW
ex_map_raw = reshape(original_maps_raw(:,1:256)',16,16,length(nan_vector));
in_map_raw = reshape(original_maps_raw(:,257:512)',16,16,length(nan_vector));
%Calculate simple difference between maps
diff_map=ex_map-in_map;
%% %% Calculate overlap of the maps ex and inh
%for all the cells
for cells = 1:length(nan_vector(incl_idx:end));
    %get the binarized maps
    bin_exc = str(nan_vector(cells)).subpixel_excMap(3:16,:);
    bin_excm= [zeros(2,16); bin_exc];
    bin_inh = str(nan_vector(cells)).subpixel_inhMap;
    %get the overlap map
    ov_map(:,:,cells) = squeeze(sum(cat(3,bin_excm,bin_inh),3));  
    ov_map_n(:,:,cells)=ov_map(:,:,cells)./max(max(ov_map(:,:,cells)));
    %Get the ex and in total per map
    ex_tot(:,cells)=str(nan_vector(cells)).excinhTotal(1);
    in_tot(:,cells)=str(nan_vector(cells)).excinhTotal(2);
end
%% Binary overlap map for EX and IN
for cells = 1:length(nan_vector(incl_idx:end));
    %get the binarized maps
    bin_exc = str(nan_vector(cells)).subpixel_excMap(3:16,:);
    bin_excm= [zeros(2,16); bin_exc]<0;
    bin_inh = str(nan_vector(cells)).subpixel_inhMap>0;
    %get the layers
   % lyr = invitro_struct(cells).layers;
    %get the overlap map
    ov_map_bin(:,:,cells) = squeeze(sum(cat(3,bin_excm,bin_inh),3))./2;  
    %layers(:,:,cells)=str(nan_vector(i)).layers;
end
% %% Morphology and Cell ID and iviv cell ID
% %morphtraces are always there but str.morph are the ones that are decent
% %traced and used for analysis
% for i=1:length(nan_vector)
%     cellID_str(i)=str(nan_vector(i)).cellID;
%     if ~isnan(str(nan_vector(i)).morph)==1;
%         morph_cells(i)=1;
%         %m_flip_a(i)=str(nan_vector(i)).morph_flip_again;
%         morph_parameters(i,:)=str(nan_vector(i)).morph;
%     else
%         morph_cells(i)=0;
%        % m_flip_a(i)=NaN;
%          morph_parameters(i,:)=ones(1,24)*NaN;
%     end
% end
% morph_cells_id=find(morph_cells==1);
% iviv_celid=find([str(nan_vector).iviv]==1);
% %% Get the morphology density maps
% for i=1:length(nan_vector)   
%         morpho_basal{i,:}=str(nan_vector(i)).morphoMap_basal_aligned;
%         morpho_apical{i,:}=str(nan_vector(i)).morphoMap_apical_aligned;     
% end
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
clean_ov=reshape(ov_map,256,cells);
clean_diff=reshape(diff_map,256,cells);
% get the number of actually populated bins
bin_num = max(pia_bins);
% allocate memory for the aligned maps
aligned_maps_ex = zeros(cell_num,16+bin_num,16);
aligned_maps_in = zeros(cell_num,16+bin_num,16);
aligned_maps_exraw = zeros(cell_num,16+bin_num,16);
aligned_maps_inraw = zeros(cell_num,16+bin_num,16);
aligned_maps_ov = zeros(cell_num,16+bin_num,16);
aligned_maps_diff = zeros(cell_num,16+bin_num,16);
% aligned_maps_basal = zeros(cell_num,16+bin_num,16);
% aligned_maps_apical = zeros(cell_num,16+bin_num,16);

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
%     % if the morpho is missing, skip
%     if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
%         aligned_maps_basal(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
%             morpho_basal{cells};        
%     end
%     % if the morpho is missing, skip
%     if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
%         aligned_maps_apical(cells,(1:16)+(1+bin_num-cell_soma),:) = ...
%             morpho_apical{cells};        
%     end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_ov = reshape(aligned_maps_ov,cell_num,[]);
aligned_maps_diff = reshape(aligned_maps_diff,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
% aligned_maps_basal = reshape(aligned_maps_basal,cell_num,[]);
% aligned_maps_apical = reshape(aligned_maps_apical,cell_num,[]);
%% Align the exc and inh maps horizontally
clean_ex = aligned_maps_ex';
clean_in = aligned_maps_in';
clean_exraw = aligned_maps_exraw';
clean_inraw = aligned_maps_inraw';
clean_ov = aligned_maps_ov';
clean_diff = aligned_maps_diff';
% clean_basal = aligned_maps_basal;
% clean_apical = aligned_maps_apical;
% get the center of the soma

soma_centers = cat(1,str(cell_idx).somaCenter);
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
    %hbin_num=6;
    % allocate memory for the aligned maps
    aligned_maps_ex = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_in = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_exraw = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_inraw = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_ov = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_diff = zeros(cell_num,16+bin_num, 16+hbin_num);
%     aligned_maps_basal = zeros(cell_num,16+bin_num, 16+hbin_num);
%     aligned_maps_apical = zeros(cell_num,16+bin_num, 16+hbin_num);
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
%          if ~isempty(morpho_basal{cells}) && sum(morpho_basal{cells}(:) > 0)
%             aligned_maps_basal(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_basal(cells,:,:);
%          end
% %         % if the morpho is missing, skip
%         if ~isempty(morpho_apical{cells}) && sum(morpho_apical{cells}(:) > 0)
%             aligned_maps_apical(cells, :, (1:16)+(1+hbin_num-cell_soma)) = clean_apical(cells,:,:);
%          end
    end
 end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
aligned_maps_ov = reshape(aligned_maps_ov,cell_num,[]);
aligned_maps_diff = reshape(aligned_maps_diff,cell_num,[]);
% aligned_maps_basal = reshape(aligned_maps_basal,cell_num,[]);
% aligned_maps_apical = reshape(aligned_maps_apical,cell_num,[]);
pia_input=pia_input(incl_idx:end);
%% %% Run 2 separate PCAs on ex and in and basal and apical morpho
%which maps to include: 
%incl_idx=65;
[coeff_ex,score_ex,latent_ex,~,explained_ex,mu] = pca(aligned_maps_ex(incl_idx:end,:,:));
[coeff_in,score_in,latent_in,~,explained_in,mu] = pca(aligned_maps_in(incl_idx:end,:,:));
%  [coeff_ex,score_ex,latent_ex,~,explained_ex,mu] = pca(aligned_maps_ex(find(slice_ori==1),:,:));
%  [coeff_in,score_in,latent_in,~,explained_in,mu] = pca(aligned_maps_in(find(slice_ori==1),:,:));
[coeff_exraw,score_exraw,latent_exraw,~,explained_exraw,mu] = pca(aligned_maps_exraw(incl_idx:end,:,:));
[coeff_inraw,score_inraw,latent_inraw,~,explained_inraw,mu] = pca(aligned_maps_inraw(incl_idx:end,:,:));
[coeff_ov,score_ov,latent_ov,~,explained_ov,mu] = pca(aligned_maps_ov(incl_idx:end,:,:));
[coeff_diff,score_diff,latent_diff,~,explained_diff,mu] = pca(aligned_maps_diff(incl_idx:end,:,:));
[coeff_com,score_com,latent_com,~,explained_com,mu] = pca([cat(1,aligned_maps_ex) cat(1,aligned_maps_in)]);
% [coeff_basal,score_basal,latent_basal,~,explained_basal,mu] = pca(aligned_maps_basal(morph_cells_id,:,:));
% [coeff_apical,score_apical,latent_apical,~,explained_apical,mu] = pca(aligned_maps_apical(morph_cells_id,:,:));
%% 
%[coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(morph_parameters(morph_cells_id,1:21));
%% Display variance explained for ex and in maps
var_exp(explained_in,explained_ex,{'Inhibition','Excitation'});legend boxoff;
%% %% Display variance explained combined
var_exp(explained_com,[],[]); 
%% Display coefficent of PCs ALIGNED in and ex
coeff_display(coeff_ex,coeff_in,bin_num,hbin_num);
%% Display coefficent of PCs ALIGNED combined
coeff_display(coeff_com(1:352,:),coeff_com(353:end,:),bin_num,hbin_num);
%% 

com=[score_ex(:,1:3) score_in(:,1:3) score_com(:,1:3)]; 
G=correlation_matrix(com,0);

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 150]);imagesc(G(7:9,1:6));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:9]);yticks([1:1:9]);
xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}'});xtickangle(45);
yticklabels({'PC1_{com}','PC2_{com}','PC3_{com}'});
c=colorbar;c.Label.String = 'R';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)

%% Plot scores value
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 600])
subplot(2,1,1)
[hN hP]=barExtended_SW(gca,score_ex(:,1:3));
 hN(1).FaceColor=[0.3 0.3 0.3];hN(2).FaceColor=[1 1 1];hN(3).FaceColor=[0.6 0.6 0.6];
  hP(1).FaceColor=[0.3 0.3 0.3];hP(2).FaceColor=[1 1 1];hP(3).FaceColor=[0.6 0.6 0.6];
 xlabel('Cells');box off;ylabel('PC score');xlim([0 148]);xticks([1:9:147])
 title('EX','Color','r');legend('PC1','PC2','PC3')
 subplot(2,1,2)
[hN hP]=barExtended_SW(gca,score_in(:,1:3));
 hN(1).FaceColor=[0.3 0.3 0.3];hN(2).FaceColor=[1 1 1];hN(3).FaceColor=[0.6 0.6 0.6];
  hP(1).FaceColor=[0.3 0.3 0.3];hP(2).FaceColor=[1 1 1];hP(3).FaceColor=[0.6 0.6 0.6];
 xlabel('Cells');box off;ylabel('PC score');xlim([0 148]);xticks([1:9:147])
 title('IN','Color','b');


end