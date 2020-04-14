%Load structure 
str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% Define array names and get the arrays without NaNs etc.
field_list = {'subpixel_excMap','subpixel_inhMap','pialD'};
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
        if contains(field_list{fields}, 'subpixel_excMap')
            field_cell{fields} = field_cell{fields}./min(field_cell{fields});
        elseif contains(field_list{fields}, 'subpixel_inhMap')
            field_cell{fields} = field_cell{fields}./max(field_cell{fields});
        end       
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});    
    
    for fields = 1:field_number
        field_cell_raw{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'subpixel_excMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        elseif contains(field_list{fields}, 'subpixel_inhMap')
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
%% Morphology and Cell ID and iviv cell ID
%morphtraces are always there but str.morph are the ones that are decent
%traced and used for analysis
for i=1:length(nan_vector)
    cellID_str(i)=str(nan_vector(i)).cellID;
    if ~isnan(str(nan_vector(i)).morph)==1;
        morph_cells(i)=1;
        %m_flip_a(i)=str(nan_vector(i)).morph_flip_again;
        morph_parameters(i,:)=str(nan_vector(i)).morph;
    else
        morph_cells(i)=0;
       % m_flip_a(i)=NaN;
         morph_parameters(i,:)=ones(1,24)*NaN;
    end
end
morph_cells_id=find(morph_cells==1);
iviv_celid=find([str(nan_vector).iviv]==1);
%% Get the morphology density maps
for i=1:length(nan_vector)   
        morpho_basal{i,:}=str(nan_vector(i)).morphoMap_basal_aligned;
        morpho_apical{i,:}=str(nan_vector(i)).morphoMap_apical_aligned;     
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
clean_exraw = aligned_maps_exraw';
clean_inraw = aligned_maps_inraw';
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
    %hbin_num=6;
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
%% Histogram of pia distribution with color coded in vivo morpho by itself
figure;set(gcf,'color','w');
histogram(pia_input,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
hold on;histogram(pia_input(morph_cells_id),'FaceColor','m','Orientation','horizontal');
hold on;histogram(pia_input(iviv_celid),'FaceColor','g','Orientation','horizontal');
legend([' Input (n=' num2str(length(pia_input)),')'],[' Morph. (n=' num2str(length(pia_input(morph_cells_id))),')']...
    ,[' In vivo (n=' num2str(length(pia_input(iviv_celid))),')']);
legend boxoff  
set(gca,'Ydir','reverse');yticks([100:100:400]);
%histogram(pia_input,'FaceColor','k','FaceAlpha',0.1)
set(gca,'FontSize',14);
%set(gca,'FontWeight','bold')
%% Setup A and Setup B
setups=[zeros(47,1);ones(100,1)];
%% Slice orientation
slice_ori=[str(nan_vector).sliceOri];
%% Plot average maps all
plot_avg_maps(str,nan_vector,ex_map,in_map,pia_input,10,0,[]);
%% 
 % plot_maps(str,ex_map_raw,in_map_raw,nan_vector,80,pia_input);
plot_maps(str,ex_map_raw,in_map_raw,nan_vector,99,pia_input);

%% Run 2 separate PCAs on ex and in and basal and apical morpho
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
[coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(morph_parameters(morph_cells_id,1:21));
%% Display variance explained for ex and in maps
var_exp(explained_in,explained_ex,{'Inhibition','Excitation'});legend boxoff  
%% %% Display variance explained combined
var_exp(explained_com,[],[]); 
%% %% Display variance explained diff maps
var_exp(explained_diff,[],[]); 
%% %% Display variance explained basal and apical
var_exp(explained_basal,explained_apical,{'Basal','Apical'}); 
%% %% Display variance explained morph
var_exp(explained_morph,[],[]); 

%COEFFCIENTS
%% Display coefficent of PCs ALIGNED in and ex
coeff_display(coeff_ex,coeff_in,bin_num,hbin_num);
%% Display coefficent of PCs ALIGNED combined
coeff_display(coeff_com(1:352,:),[],bin_num,hbin_num);
%% Display coefficent of PCs ALIGNED diff
coeff_display(coeff_diff,[],bin_num,hbin_num);
%% Display correlation between all PCs Multiple comparison 
com=[score_ex(:,1:3) score_in(:,1:3) score_com(:,1:3)]; 
correlation_matrix(com,1);
xticks([1:1:9]);yticks([1:1:9]);
xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','PC1_{com}','PC2_{com}','PC3_{com}'});xtickangle(45)
yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','PC1_{com}','PC2_{com}','PC3_{com}'});ytickangle(45);xlabel('Feature');ylabel('Feature');
%% %%Display desired correlations between PCs/pia
corr_plot(score_ex(:,1),score_com(:,1),[],{'PC1ex','PC1com','Pial depth'});
%% Load fraction and absolute input of maps for ex and inh
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv L23h L4h L5h pialD layer_assign] = iviv_profiles(nan_vector(incl_idx:end),str);
frac_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv];
abs_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv];
layer_assign=[zeros(length(nan_vector(incl_idx:end)),2) layer_assign'];
%% Normalizied maps to calculate abs
[frac_exh_n abs_exh_n frac_inh_n abs_inh_n frac_exv_n abs_exv_n frac_inv_n abs_inv_n] = iviv_profiles_norm(nan_vector,ex_map,in_map);
frac_exv_n=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv_n];
abs_exv_n=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv_n];
%% Raw Max per map EX and IN
for i=1:length(nan_vector)
max_map_ex(i)=abs([min(min(ex_map_raw(:,:,i)))]);
max_map_in(i)=abs([max(max(in_map_raw(:,:,i)))]);
end
%% Calculate the fraction from the overlap for verti and hori
[frac_ovh abs_ovh frac_ovv abs_ovv] = iviv_profiles_ov(ov_map);
%% Calculate fraction using the overlap binary maps
[frac_bh abs_bh frac_bv abs_bv] = iviv_profiles_ov(ov_map_bin);
%% Calculate difference between ex and in (ex-in) USING THIS FOR NOW
frac_diffv=frac_exv_m-frac_inv;
frac_diffh=frac_exh-frac_inh;
%% Calculate VERTICAL ex and in fraction for L23, L4,  L5
L23fr=[nanmean(frac_exv_m(:,3:5),2) nanmean(frac_inv(:,3:5),2)];
L4fr=[nanmean(frac_exv_m(:,6:7),2) nanmean(frac_inv(:,6:7),2)];
L5fr=[nanmean(frac_exv_m(:,8:10),2) nanmean(frac_inv(:,9:11),2)];
L23frov=[nanmean(frac_ovv(:,3:5),2)];
L4frov=[nanmean(frac_ovv(:,6:7),2)];
L5frov=[nanmean(frac_ovv(:,8:10),2)];
%% Calculae difference between fraction 
L5frt=L5fr;
L5frt(find(L5fr(:,1)==0),1)=NaN ;
L5frt(find(L5fr(:,2)==0),2)=NaN ;
diffL23fr=L23fr(:,1)-L23fr(:,2);
diffL4fr=L4fr(:,1)-L4fr(:,2);
diffL5fr=L5frt(:,1)-L5frt(:,2);
%% Calculate fraction using the layer assignement
for i=1:length(frac_exv_m)
    L4fr_l(i,:)=[nanmean(frac_exv_m(i,find(layer_assign(i,:)==3))) nanmean(frac_inv(i,find(layer_assign(i,:)==3)))];   
end
%% Display alignement startegy for layers
figure;set(gcf,'color','w');imagesc(layer_assign')
hold on;line([47 47], [1 16],'Color','k','LineStyle','--');set(gca,'FontSize',10);
ylim([1 16]);yticks([1:1:16]);ylabel('Rows');xlabel('Cells');%c=colorbar;
%% Calculate difference between ex and in (ex-in) medial and lateral for whole map
%mean
frh_medial=[nanmean(frac_exh(:,1:8),2) nanmean(frac_inh(:,1:8),2)];
frh_lateral=[nanmean(frac_exh(:,9:end),2) nanmean(frac_inh(:,9:end),2)];
frh_diff_medial=nanmean(frac_diffh(:,1:8),2);
frh_diff_lateral=nanmean(frac_diffh(:,9:end),2);
%sum
frh_medial_s=[sum(frac_exh(:,1:8),2) sum(frac_inh(:,1:8),2)];
frh_lateral_s=[sum(frac_exh(:,9:end),2) sum(frac_inh(:,9:end),2)];
%% Calculate the maximum horizontal span overall
for i=1:length(nan_vector(incl_idx:end))
tmp=find(frac_exh(i,:)>0);
ex_spanh(i)=tmp(end)-tmp(1);
tmp=[];
tmp=find(frac_inh(i,:)>0);
in_spanh(i)=tmp(end)-tmp(1);
tmp=[];
end
%% Calculate the maximum horizontal span overall per layer
[ex_spanhL23 ex_spanhL4 ex_spanhL5 in_spanhL23 in_spanhL4 in_spanhL5] = span_perLayer(ex_map,in_map,nan_vector);
%% total sum absolute for ex and inh in pA
for i=1:length(nan_vector(incl_idx:end))
%whole map
tmp1=ex_map(3:end,:,i);
tmp2=in_map(:,:,i);
% tmp1=ex_map(:,:,i);
% tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
%tot_input(i,:)=[sum(tmp1(:))/length(tmp1);sum(tmp2(:))/length(tmp2)];
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map(3:5,:,i);
tmp2=in_map(3:5,:,i);
%  tmp1=ex_map(3:5,:,i);
%  tmp2=in_map(3:5,:,i);
tot_inputL23(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map(6:7,:,i);
tmp2=in_map(6:7,:,i);
%  tmp1=ex_map(6:7,:,i);
%  tmp2=in_map(6:7,:,i);
tot_inputL4(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
end
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,[]);
%% Display fraction for ex and in as well as overlap for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_ovv*-1,frac_ovh*-1,[]);
%% Display fraction for ex and in as well as overlap bin for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_bv*-1,frac_bh*-1,[]);
%% Display both maximum horizontal span and diff betwen ex and in COLUMN APPROACH
xf=69;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 450, 500]);
subplot(3,2,1);
;hold on;
h4 = histogram(diffL23fr,8);h4.BinWidth = 0.025;h4.EdgeColor = 'k';h4.FaceColor = 'm';hold on;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylim([0 50]);yticks([0:25:50])
title('Vertical');ylabel('Cell counts');hold on;set(gca,'FontSize',10);
xlim([-0.3 0.3]);text(-0.2,50,'IN','Color','b');text(0.2,50,'EX','Color','r');text(-0.2,40,'L2/3','Color','k');
hold on;plot(nanmean(diffL23fr),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');
subplot(3,2,3);
h4 = histogram(diffL4fr,8);h4.BinWidth =  0.025;h4.EdgeColor = 'k';h4.FaceColor = 'g'
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylabel('Cell counts');hold on;set(gca,'FontSize',10);box off;
text(-0.2,40,'L4');
hold on;plot(nanmean(diffL4fr),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
xlim([-0.3 0.3]);
subplot(3,2,5);
h4 = histogram(nonzeros(diffL5fr),8);h4.BinWidth =  0.025;h4.FaceColor = [0.5 0.5 0.5];ylabel('Cell counts');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');
xlim([-0.3 0.3]);xlabel('\Delta EX-IN vertical fraction');hold on;set(gca,'FontSize',10);
text(-0.2,40,'L5','Color','k');
hold on;plot(nanmean(diffL5fr),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
%Horizontal span
subplot(3,2,2);
h4 = histogram((ex_spanhL23-in_spanhL23)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'm';
xlim([-10*xf 10*xf]);hold on; title('Horizontal');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,50,'IN','Color','b');text(7*xf,50,'EX','Color','r');
text(-7*xf,40,'L2/3','Color','k');
hold on;plot(nanmean((ex_spanhL23-in_spanhL23)*xf),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,4);
hold on;h4 = histogram((ex_spanhL4-in_spanhL4)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'g';
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,40,'L4','Color','k');
hold on;plot(nanmean((ex_spanhL4-in_spanhL4)*xf),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,6);
hold on;h4 = histogram((ex_spanhL5-in_spanhL5)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];;
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
xlabel('\Delta EX-IN horizontal span (µm)');hold on;set(gca,'FontSize',10);
text(-7*xf,40,'L5','Color','k');
hold on;plot(nanmean((ex_spanhL5-in_spanhL5)*xf),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
%group based on pial depth
g1=find(pia_input>=220);g2=find(pia_input<220);g3=[];
gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
gv(g1)=1;gv(g2)=2;gv(g3)=3;
%call function
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
g1=[];g2=[];g3=[];
%% Pial depth correlation plot with fraction 
plot_fraction(L23fr,L4fr,L5fr,pia_input);
%% Get the somax and y coordinates using the grid pattern
for i=1:length(nan_vector)
    if slice_ori(i)==1
somax1(i)=str(nan_vector(i)).somaCenter(1);
    else slice_ori(i)==0
    somax1(i)=-str(nan_vector(i)).somaCenter(1);    
    end
somay1(i)=552-str(nan_vector(i)).somaCenter(2);
end

for i=1:length(nan_vector)
   
somax(i)=str(nan_vector(i)).subpixel_soma(1);    
  
somay(i)=552-str(nan_vector(i)).subpixel_soma(2);
end
%% Calculate centroid of map and angle towards it 
%fake map
somax=somax1
map_fake=zeros(16,16,cells)
map_fake(4,9,:)=1;
%IN ang centroid
[out_ang_in] = centroid_map(in_map(:,:,:),somax,pia_input,[1:cells],0);
[out_ang_inL23] = centroid_map(in_map(3:5,:,:),somax,pia_input,[1:cells],2);
[out_ang_inL4] = centroid_map(in_map(6:7,:,:),somax,pia_input,[1:cells],5);
[out_ang_inL5] = centroid_map(in_map(8:10,:,:),somax,pia_input,[1:cells],7);
%EX ang centroid
[out_ang_ex] = centroid_map(ex_map(:,:,:),somax,pia_input,[1:cells],0);
[out_ang_exL23] = centroid_map(ex_map(3:5,:,:),somax,pia_input,[1:cells],2);
[out_ang_exL4] = centroid_map(ex_map(6:7,:,:),somax,pia_input,[1:cells],5);
[out_ang_exL5] = centroid_map(ex_map(8:10,:,:),somax,pia_input,[1:cells],7);
%difference map
[out_ang_diff] = centroid_map(diff_map(3:5,:,:),somax,pia_input,[1:cells],0);
%fake map
[out_ang_fake] = centroid_map(map_fake(3:5,:,:),somax,pia_input,[1:cells],2);
%% Plotting the centroid with vector pointing towards it
ang1=out_ang_exL23;
ang2=out_ang_inL23;
ang3=out_ang_exL4;
idxtp=1:147;
%call function quiver
quiver_centroid(ang1,ang2,ang3,idxtp,ex_map,in_map,128,2);
%% Distributions of actual centroid all input cells for EX and IN
a=1:147;
centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,0,[],{});
%% Histogram centroid X
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 600, 250]);
subplot(1,3,1);
h1=histogram(out_ang_exL23(:,3)*69-out_ang_inL23(:,3)*69)
h1.EdgeColor = 'k';h1.FaceColor = 'm';hold on;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');
box off;
xlim([-200 200]);yticks([0:20:60])
ylabel('Counts');
hold on;title('L23');set(gca,'FontSize',10);
subplot(1,3,2);
h1=histogram(out_ang_exL4(:,3)*69-out_ang_inL4(:,3)*69)
hold on;title('L4');set(gca,'FontSize',10);
h1.EdgeColor = 'k';h1.FaceColor = 'g';hold on;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');
xlim([-200 200]);;xlabel('\Delta EX-IN Horizontal')
box off;yticks([0:20:60])
subplot(1,3,3);
h1=histogram(out_ang_exL5(:,3)*69-out_ang_inL5(:,3)*69)
hold on;title('L5');set(gca,'FontSize',10);
box off;
h1.EdgeColor = 'k';h1.FaceColor = [0.5 0.5 0.5];hold on;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');
xlim([-400 400])
box off;yticks([0:20:60])
%% % Relation betwen ex and in centroid
ang1=out_ang_exL23;
ang2=out_ang_inL23;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 900, 200]);
subplot(1,3,1);
pointsize=20;
[cmap]=buildcmap('kmgy');
scatter(ang1(:,3)*69-ang2(:,3)*69,ang1(:,4)*69-ang2(:,4)*69,pointsize,pia_input,'filled');set(gcf,'color','w');xlabel('Ex-In Horizontal');ylabel('Ex-In Vertical');set(gca,'FontSize',10)
hold on;line([-200 200], [0 0],'Color','k','LineStyle','--');hold on;line([0 0], [-200 200],'Color','k','LineStyle','--'); c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (µm)';
xticks([-200:100:200]);yticks([-200:100:200]);set(gca,'FontSize',10);axis square;
subplot(1,3,2);
scatter(ang1(:,5),ang2(:,5),pointsize,pia_input,'filled');set(gcf,'color','w');;xlabel('Excitation (deg)');ylabel('Inhibition (deg)');c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (µm)';hold on; refline(1,0);set(gca,'FontSize',10);axis square;
subplot(1,3,3);
scatter(ang1(:,8)*69,ang2(:,8)*69,pointsize,pia_input,'filled');set(gcf,'color','w');xlabel('Excitation');ylabel('Inhibition ');c=colorbar;colormap(cmap);c.Label.String = 'Pial depth (µm)';hold on;refline(1,0);set(gca,'FontSize',10);axis square;
%% Difference between EX and IN centroid
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 250, 600]);
subplot(3,1,1);
plot(out_ang_exL23(:,3)*69-out_ang_inL23(:,3)*69,out_ang_exL23(:,4)*69-out_ang_inL23(:,4)*69,'o','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on;line([-200 200], [0 0],'Color','k','LineStyle','--');hold on;line([0 0], [-200 200],'Color','k','LineStyle','--');box off
;ylabel('\Delta EX-IN Vertical')
hold on;title('L23');set(gca,'FontSize',10);
subplot(3,1,2);
plot(out_ang_exL4(:,3)*69-out_ang_inL4(:,3)*69,out_ang_exL4(:,4)*69-out_ang_inL4(:,4)*69,'o','MarkerFaceColor','g','MarkerEdgeColor','k')
hold on;line([-200 200], [0 0],'Color','k','LineStyle','--');hold on;line([0 0], [-200 200],'Color','k','LineStyle','--');box off
;ylabel('\Delta EX-IN Vertical')
hold on;title('L4');set(gca,'FontSize',10);
subplot(3,1,3);
plot(out_ang_exL5(:,3)*69-out_ang_inL5(:,3)*69,out_ang_exL5(:,4)*69-out_ang_inL5(:,4)*69,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k')
hold on;line([-400 400], [0 0],'Color','k','LineStyle','--');hold on;line([0 0], [-200 200],'Color','k','LineStyle','--');box off
xlabel('\Delta EX-IN Horizontal');ylabel('\Delta EX-IN Vertical');xlim([-400 400])
hold on;title('L5');
set(gca,'FontSize',10);
%% Ex vs IN with refline
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 900, 200]);
subplot(1,3,1);
plot(out_ang_exL23(:,3)*69-out_ang_exL23(:,1)*69,out_ang_inL23(:,3)*69-out_ang_inL23(:,1)*69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);
ref.Color='k';box off;xlim([-200 200]);ylim([-200 200]);hold on;title('L23');xticks([-200:100:200]);yticks([-200:100:200]);
subplot(1,3,2);
plot(out_ang_exL4(:,3)*69-out_ang_exL4(:,1)*69,out_ang_inL4(:,3)*69-out_ang_inL4(:,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);
ref.Color='k';box off;xlim([-200 200]);ylim([-200 200]);hold on;title('L4');xticks([-200:100:200]);yticks([-200:100:200]);
subplot(1,3,3);
plot(out_ang_exL5(:,3)*69-out_ang_exL5(:,1)*69,out_ang_inL5(:,3)*69-out_ang_inL5(:,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);;
ref.Color='k';box off;xlim([-400 400]);ylim([-400 400]);hold on;title('L5');xticks([-400:200:400]);yticks([-400:200:400]);
%% Test
par1=out_ang_exL23(:,3)*69-out_ang_exL23(:,1)*69;
par2=out_ang_inL23(:,3)*69-out_ang_inL23(:,1)*69;
find(~isnan(par2))
P = polyfit(par1(find(~isnan(par2))),par2(find(~isnan(par2))),1);
    yfit = P(1)*par1(find(~isnan(par2)))+P(2);
 yfit2=1*par1(find(~isnan(par2)));
 
 [r t]=ttest2(yfit,yfit2)
%% Desired correlation between angle and inhbition PC2in
corr_plot(abs(out_ang_inL23(:,5)),score_ex(:,3),pia_input,{'Angle centroid L2/3 IN','PC2_{in}','Pial depth (µm)'});ylabel('PC2_{in}','Color','b');xlabel('Angle centroid L2/3 IN','Color','b');
%% Desired correlation between horizontal centroid and PC3in
corr_plot((out_ang_inL23(:,3)-out_ang_inL23(:,1))*69,score_in(:,3),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});xlim([-150 150]);ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% Fraction medial/lateral with scores PCin3
corr_plot(frh_medial_s(:,2),score_in(:,3),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% Fraction medial/lateral with centroid X IN
corr_plot(frh_medial_s(:,2),(out_ang_inL23(:,3)-out_ang_inL23(:,1))*69,pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% 
corr_plot((out_ang_inL23(:,3)-out_ang_inL23(:,1))*69,score_in(:,3),[],{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});ylabel('PC3_{in}','Color','b');xlabel('Horizontal centre of mass L2/3','Color','b');set(gca,'FontSize',12)


%% 
corr_plot(in_spanhL4',score_in(:,2),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'})
%% Pial depth PC3ex
corr_plot(score_ex(:,1),pia_input,[],{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});xlabel('PC1_{ex}','Color','r');ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',12);set(gca,'Ydir','reverse');
%% 
corr_plot(90-abs((out_ang_inL23(:,5))),pia_input,[],{'Horizontal centroid L2/3 IN','PC3_{in}','Pial depth'});ylabel('Pial depth (µm)','Color','k');xlabel('Angle Centre of mass L23','Color','b');set(gca,'FontSize',12);set(gca,'Ydir','reverse');

%% FAKE MAPS CENTROID
corr_plot(out_ang_fake(:,5),score_in(:,2),pia_input,{'Angle centroid L23 fake','PCin2','Pial depth'});set(gca,'Ydir','reverse');
%% Display correlation between all PCs Multiple comparison 
com=[score_ex(:,1:3) score_in(:,1:3) L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) out_ang_exL23(:,5)  out_ang_exL4(:,5) out_ang_inL23(:,5)  out_ang_inL4(:,5) frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2)]; 
G=correlation_matrix(com,1);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 500]);imagesc(G(7:18,1:6));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:18]);yticks([1:1:18]);
% xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});xtickangle(45)
% yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});ylabel('Feature');
xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}'});xlabel('Feature');xtickangle(45);set(gca,'FontSize',12)
yticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});ylabel('Feature');ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% More to compare between PCs and actual data
com=[score_ex(:,1:3) score_in(:,1:3) L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) out_ang_exL23(:,5)  out_ang_exL4(:,5) out_ang_inL23(:,5)  out_ang_inL4(:,5) frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) ex_spanhL23' in_spanhL23' ex_spanhL4' in_spanhL4' ex_spanhL5' in_spanhL5' ex_spanh' tot_inputL4 sum(abs_exv_n(:,6:7),2) ]; 
G=correlation_matrix(com,1);

%% PC2ex with abs exv_n
corr_plot(sum(abs_exv_m(:,6:7),2),score_ex(:,2),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});ylabel('PC2_{ex}','Color','r');xlabel('Sum normalized input L4');
%% 


%MORPHOLOGY AND INPUT
%% Get 16x16 maps for basal and apical
for i=1:length(nan_vector)
     if isempty(morpho_basal{i,:,:})==0
     ba_map(:,:,i)=morpho_basal{i,:,:};
     ap_map(:,:,i)=morpho_apical{i,:,:};
     else
     ba_map(:,:,i)=ones(16,16)*NaN;
     ap_map(:,:,i)=ones(16,16)*NaN;
     end
end
%% Max density + total sum per cell apical and basal
max_densba=[max(max(ba_map(:,:,:)))];
max_densba=reshape(max_densba,1,147);
max_densap=max(max(ap_map(:,:,:)));
max_densap=reshape(max_densap,1,147);
%SUM of ap and bas overall
for i=1:length(nan_vector)
    sum_densap(i)=sum(sum(ap_map(:,:,i)));
    sum_densba(i)=sum(sum(ba_map(:,:,i)));
end
%% split in medial and lateral
for i=1:length(nan_vector)
MLba_diff(i)=(sum(nonzeros(ba_map(:,1:8,i)))/length(nonzeros(ba_map(:,1:8,i))))-(sum(nonzeros(ba_map(:,9:end,i)))/length(nonzeros(ba_map(:,9:end,i))));
MLap_diff(i)=(sum(nonzeros(ap_map(:,1:8,i)))/length(nonzeros(ap_map(:,1:8,i))))-(sum(nonzeros(ap_map(:,9:end,i)))/length(nonzeros(ap_map(:,9:end,i))));
end
%% Morphology parameters vs PCs
com=[];com=[score_ex(:,1:3) score_in(:,1:3) morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff' sum_densap' sum_densba']; 
correlation_matrix(com,1);title('Morphology vs PC input scores');
%% Read out map correlation with morphology density
corr_exc_apical=[str(:).corr_exc_apical];
corr_exc_basal=[str(:).corr_exc_basal];
corr_inh_apical=[str(:).corr_inh_apical];
corr_inh_basal=[str(:).corr_inh_basal];
%% %% Morphology parameters vs real parameters vertical
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) out_ang_inL23(:,5) out_ang_inL23(:,3) morph_parameters(:,1:21) sum_densap' sum_densba']; 
correlation_matrix(com,1);title('Morphology vs Vertical input');
%% Total input apical dendrite branch points
corr_plot(morph_parameters(:,4),in_spanhL23',[],{'Branch points','Total input EX','Pial depth (µm)'})
%% Total input apical branch  density
corr_plot(sum_densap',abs(tot_input(:,1)),pia_input,{'Branch points','Total input EX','Pial depth (µm)'})
%% %% Morphology parameters vs real parameters horizontal
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff' sum_densap' sum_densba']; 
correlation_matrix(com,0);title('Morphology vs Horizontal input');



%% %%%%%%%%%%%%%%%%

%% In vivo features alone

%% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 [od_out sftf_out sftf_sel sftf_pref spon_out pia_all] = concat_invivo(L23_PC);
%%  
com=[od_out(:,:) sftf_pref(:,10:18) pia_all']; 
correlation_matrix(com,0);title('In vivo alone');
 %% %%Display desired correlations between pial depth and L23/L4
corr_plot(com(:,1),abs(com(:,8)),pia_all,{'OSI bino','OSI SFTF','Pial depth'});
%% Correlation all
com=[];
com=[sftf_sel]; 
correlation_matrix(com,1);title('In vivo alone ORI');
 %% %%Display desired correlations between pial depth and L23/L4
corr_plot(com(:,1),abs(com(:,7)),com(:,10),{'L4ex input','Apical width/height','XSA'});
 %% %%Display desired correlations between pial depth and L23/L4
corr_plot(sftf_out(find(sftf_out(:,3)>0.01),5),od_out(find(sftf_out(:,3)>0.01),8),sftf_out(find(sftf_out(:,3)>0.01),2),{'L4ex input','Apical width/height','XSA'});
%% Discretize based on pia
discretize_plot(pia_all,3,od_out(:,3),1);xlabel('Pial depth bins');ylabel('Orientation preference');xticks([1:1:3]);title('all in vivo cells');ylim([0 180])
%% 
discretize_plot(pia_all(find(od_out(:,3)>=90)),3,od_out(find(od_out(:,3)>=90),3),1);
%% 
discretize_plot(pia_all,3,abs(sftf_pref(:,21)),1);
%% 
discretize_plot(pia_all,3,abs(od_out(:,3)),1);ylim([0 180]);xticklabels({''});
hold on;line([1 3], [90 90],'Color','k','LineStyle','--');
%% 
discretize_plot(pia_input,3,abs(od_out_iviv(:,4)),1);ylim([0 180]);xticklabels({''});
hold on;line([1 3], [90 90],'Color','k','LineStyle','--');



%%

%%  IN VIVO IN VITRO
%% Read out paramteres from iviv structure
[od_out_iviv spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(str,nan_vector);
%% Get IDS for ORI
id_ori=find(~isnan(od_out_iviv(:,4)));
%% %read out noise correlations
noise_corr=[str(:).noise];
%% Correlations PCs input and OD out
ia=find(od_out_iviv(:,2)>0.3)  
%b=find(od_out_iviv(:,1)>0.3)  
%[ia ib]=intersect(a,b)
com=[];com=[score_ex(ia,1:3) score_in(ia,1:3) od_out_iviv(ia,:) noise_corr(ia)']; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations PCs input and OD out COM PCs
ia=find(od_out_iviv(:,2)>0.3)  
%b=find(od_out_iviv(:,1)>0.3)  
%[ia ib]=intersect(a,b)
com=[];com=[score_com(ia,1:3) od_out_iviv(ia,:)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations PCs input and sftft out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) sftf_out_iviv(:,2:end)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations Input vertical vs ODout
%a=find(od_out_iviv(:,1)>0)  
a=find(noise_corr<0.1)
com=[];com=[L23fr(a,1) L4fr(a,1) L5fr(a,1) L23fr(a,2) L4fr(a,2) L5fr(a,1) tot_input(a,1) tot_input(a,2) out_ang_inL23(a,3)  sum(abs_exv_n(a,1:5),2) ex_spanhL23(a)' ex_spanhL4(a)' frh_lateral_s(a,1) frh_lateral_s(a,2) od_out_iviv(a,1:8) noise_corr(a)']; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Score in with Orientation
a=find(od_out_iviv(:,1)>0.3)  
b=find(pia_input<300)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(ia,4),score_in(ia,2),pia_input(ia),{'11','2','pia'})
corr_plot(od_out_iviv(ia,4),out_ang_inL23(ia,5),pia_input(ia),{'Orientation','Angle centre of mass L23 in','Pial depth'})
;set(gca,'FontSize',12);
%% Correlations Input vertical vs ODout
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) tot_input(:,1) od_out_iviv]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Span vs od out
a=find(od_out_iviv(:,2)>0.3)  
com=[];com=[ex_spanhL23(a)' in_spanhL23(a)' ex_spanhL4(a)' in_spanhL4(a)'  abs(ex_spanhL23(a)'-in_spanhL23(a)')*69 od_out_iviv(a,:)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Correlations Input vertical vs SFTF pref
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) sftf_out_pref_iviv(:,1:3)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% 
%corr_plot(od_out_iviv(a,[4]),abs(out_ang_inL23(a,5)),pia_input(a),{'L4ex input','Apical width/height','XSA'});
corr_plot(od_out_iviv(a,[4]),abs(out_ang_inL23(a,3)-out_ang_inL23(a,1))*69,pia_input(a),{'L4ex input','Apical width/height','XSA'});

%% Display correlation between all input maps and in vivo
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) abs(out_ang_exL23(:,5)) abs(out_ang_inL23(:,5)) abs(out_ang_exL23(:,3)-out_ang_exL23(:,1)) abs(out_ang_inL23(:,3)-out_ang_inL23(:,1)) od_out_iviv(:,[1 2 3 4 5 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
%% Subsample OSI DSI ODI TW Ca peak
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) abs(out_ang_exL23(:,5)) abs(out_ang_inL23(:,5)) abs(out_ang_exL23(:,3)-out_ang_exL23(:,1)) abs(out_ang_inL23(:,3)-out_ang_inL23(:,1)) pia_input od_out_iviv(:,[1 2 3 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:14]);yticks([1:1:14]);
%% Display correlation between all input maps and in vivo
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 300]);imagesc(G(10:end,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:13]);yticks([1:1:5]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','TW','Ca_{peak}'});set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% Subsample data OSI
a=find(od_out_iviv(:,1)>0.25)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) abs(out_ang_exL23(a,5)) abs(out_ang_inL23(a,5)) abs(out_ang_exL23(a,3)-out_ang_exL23(a,1)) abs(out_ang_inL23(a,3)-out_ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
% xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});xtickangle(45)
% yticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});ytickangle(45)
%% ORI with cutoff
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 200]);imagesc(G(13,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:1]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'ORI'});set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);
%% Subsample data DSI
a=find(od_out_iviv(:,2)>0.25)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) abs(out_ang_exL23(a,5)) abs(out_ang_inL23(a,5)) abs(out_ang_exL23(a,3)-out_ang_exL23(a,1)) abs(out_ang_inL23(a,3)-out_ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
%% DIR with cutoff
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 200]);imagesc(G(14,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:1]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'DIR'});set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);

%% IMPORTANT CALCULATIONS
%% Direction difference Lateral vs medial
g1=[];
g2=[];
a=find(od_out_iviv(:,2)>0.25)  
b=find(pia_input<300)
[ia ib]=intersect(a,b)
%par=frh_lateral_s(ia,2)-frh_medial_s(ia,2)
%par=out_ang_inL23(ia,3)
%par=L4h(ia,4)-L4h(ia,3)
par=(out_ang_in(ia,3)-out_ang_in(ia,1))*69
g1=find(od_out_iviv(ia,5)>0 & od_out_iviv(ia,5)<225) 
g2=find(od_out_iviv(ia,5)>240 & od_out_iviv(ia,5)<360) 
% tmp=zeros(length(ia),1)
% tmp(g1)=g1;
% g2=find(tmp==0);
%g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<360) 
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2])
xticklabels({'~135°','~315'})
xlabel('Direction');ylabel('Difference IN Lateral-Medial');  set(gca,'FontSize',12);
%xlabel('Direction');ylabel('\Delta IN Horizontal Centroid');  set(gca,'FontSize',12);
%% 
figure;set(gcf,'color','w');
histogram(od_out_iviv(ia,5),12);box off
ylabel('Counts')
xlabel('Direction preference')

%% Correlaion of the above
a=find(od_out_iviv(:,2)>0.3)  
b=find(pia_input<300)
%b=find(od_out_iviv(:,1)>0.3)  
[ia ib]=intersect(a,b)  
corr_plot(frh_lateral_s(ia,2)-frh_medial_s(ia,2),od_out_iviv(ia,5),od_out_iviv(ia,2),{'OSI','IN INPUT','pia'});
%% SF vs abs difference between L23 EX/IN span
a=find(abs(ex_spanhL23(:)-in_spanhL23(:))*69<200)
par=abs(ex_spanhL23(a)-in_spanhL23(a))*69;
%par=(in_spanhL23(a))*69;
%par=score_ex(a,2);
%par=in_spanhL23(a)*69;
 g1=find(sftf_out_iviv(a,2)==0.02)
 g2=find(sftf_out_iviv(a,2)==0.08)
%  g1=find(od_out_iviv(:,1)<0.5)
%  g2=find(od_out_iviv(:,1)>0.5)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'0.02','0.08'})
xlabel('SF');ylabel('| \Delta EX-IN horizontal span (µm) |');
%ylim([0 200]);
%% Scatter plot for SF finding
par1=ex_spanhL23(a)*69;
par2=in_spanhL23(a)*69;
figure;scatter(par1(g1),par2(g1),'mo','filled');set(gcf,'color','w')
hold on;scatter(par1(g2),par2(g2),'go','filled');
 refline(1,0);
 xlim([4*9 16*69]);
 ylim([4*69 16*69]);
 ylabel('IN span horizontal (µm)');
 xlabel('EX span horizontal (µm)');
  set(gca,'FontSize',12);
  legend('0.02','0.08')
%% L4 fraction with ORI pref
a=find(od_out_iviv(:,1)>0.25)
corr_plot(od_out_iviv(a,4),L4fr(a,1),[],{'Orientation preference','L4ex fraction'});ylabel('EX L4 fraction','Color','r');set(gca,'FontSize',12);
%% Ori pref with pial depth
corr_plot(od_out_iviv(a,4),pia_input(a),[],{'ORI','Pial depth'});set(gca,'FontSize',12);
%% L4 fraction with pial dpeth
corr_plot(L4fr(a,1),pia_input(a),od_out_iviv(a,4),{'L4exfraction','Pial depth','ORI'});set(gca,'FontSize',12);
%% 
corr_plot(L4fr(a,1),abs(out_ang_inL23(a,3)-out_ang_inL23(a,1))*69,pia_input(a),{'L4exfraction','Pial depth','ORI'});set(gca,'FontSize',12);
%% Plot x with L4 fraction
%corr_plot(abs(out_ang_inL23(ia),3)-out_ang_inL23(ia,1))*69,L4fr(ia,1),[],{'|\Delta Centroid L23 X|','L4exfraction','ORI'});set(gca,'FontSize',12);
%% 
corr_plot(90-abs(out_ang_inL23(:,5)),pia_input,[],{'|Angle Centroid L23|','Pial depth','ORI'});set(gca,'FontSize',12);
%% categorize cells basen on orientations
g1=[];
g2=[];
a=find(od_out_iviv(:,1)>0.25)  
b=find(pia_input<300)
[ia ib]=intersect(a,b)
%par=frh_lateral_s(ia,2)-frh_medial_s(ia,2)
%par=abs(out_ang_exL23(ia,5))
par=L23fr(ia,1)
%par=abs(out_ang_inL23(ia,3)-out_ang_inL23(ia,1))*69
%par=abs(out_ang_exL23(ia,10))*69
%par=abs(out_ang_inL23(ia,1))*69
%par=in_spanhL23(ia)*69
%par=score_in(ia,2)
%par=pia_input(ia)
g1=find(od_out_iviv(ia,4)>0 & od_out_iviv(ia,4)<90) 
% tmp=zeros(length(ia),1)
% tmp(g1)=g1;
g2=find(od_out_iviv(ia,4)>90 & od_out_iviv(ia,4)<180) 
%g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<360) 
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2])
xticklabels({'0-90°','90-180°'})
xlabel('Orientation');
%ylabel('|Angle Centroid L23|','Color','r');  set(gca,'FontSize',12);
%ylabel('|\Delta Centroid L23 X|','Color','r');  set(gca,'FontSize',12);
%ylabel('EX L4 fraction','Color','r');
%ylabel('IN L4 fraction','Color','b');
ylabel('EX L2/3 fraction','Color','r');
%ylabel('Pial depth (µm)')
set(gca,'FontSize',12);
%% 
a=find(od_out_iviv(:,1)>0.25)  
%fe=od_out_iviv(a,4)
fe=pia_input(a)
centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,1,fe,{'ORI'});
%% 
pointsize=40;
figure;
h4 = scatter(out_ang_inL23(a,3)*69-ang2(a,1)*69,out_ang_inL23(a,4)*69-out_ang_inL23(a,2)*69,pointsize,fe,'filled');xlim([-4*69 4*69]);ylim([-4*69 8*69]);
[cmap]=buildcmap('ybk');
colormap(cmap);%text(-250,50,'L2/3');
ang2=out_ang_exL4;
yticks([-200:200:600])
set(gca,'Ydir','reverse');hold on;line([0 0], [-4*69 8*69],'Color','k','LineStyle','--');hold on;line([-4*69 8*69],[0 0],'Color','k','LineStyle','--');hold on;plot(0,0,'^r');
hold on;title('IN Centre of mass');c=colorbar;c.Label.String='ORI';  
ylim([-100 100])
xlim([-100 100])
%% 

figure;scatter(par(g2)*-1,ones(length(g2),1),'og','filled')
hold on;scatter(zeros(length(g1),1),par(g1)*-1,'om','filled')
hold on;
xlim([-80 50]);
ylim([-80 50]);
%% Correlations Input horizontal vs ODout
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' od_out_iviv]; 
correlation_matrix(com,0);title('Input horizontal ODout iviv');
%% Correlations PCs input and spon out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) spon_out_iviv]; 
correlation_matrix(com,1);title('PC inputs spon_out iviv');
%% Correlations Vertical input and spon out
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) L23frov L4frov spon_out_iviv]; 
correlation_matrix(com,1);title('Input vertical spon_out iviv');
%% Correlations Horizontal input and spon out
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)'  spon_out_iviv]; 
%% Correlations PCs input and sftf_out_iviv
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov sftf_out_iviv(:,[1:6])]; 
correlation_matrix(com,0);title('PC inputs sftf_out iviv');
%% Correlations PCs input and sftf_out_iviv
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' sftf_out_iviv]; 
correlation_matrix(com,0);title('PC inputs sftf_out iviv');

%% %% Intersectional approach ORI
a=find(L4fr(:,1)<0.3)  
b=find(od_out_iviv(:,1)>0.3)
%c=find(od_out_iviv(:,4)<180)
c=find(pia_input<260)
[ia ib]=intersect(a,b)
[iaa ic]=intersect(ia,c);
corr_plot(od_out_iviv(iaa,4),out_ang_inL23(iaa,5),pia_input(iaa),{'ORI','L4ex fraction','pia'});set(gca,'FontSize',12);
%% responsive vs non responsive
g1=find(sftf_out_iviv(:,1)==0);
g2=find(od_out(:,1)==1);
%c_resp=rf(find(iv_Ca(find(iv_resp2==1))>200))
gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
gv(g1)=1;gv(g2)=2;gv(g3)=3;
gv(find(gv==0))=NaN;
%call function
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
g1=[];g2=[];g3=[];

%% 

%% Ex vs IN with refline subsampled
a=find(od_out_iviv(:,1)>0.25) 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 900, 200]);
subplot(1,3,1);
plot(abs(out_ang_exL23(a,3)*69-out_ang_exL23(a,1)*69),abs(out_ang_inL23(a,3)*69-out_ang_inL23(a,1)*69),'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);
ref.Color='k';box off;xlim([-200 200]);ylim([-200 200]);hold on;title('L23');xticks([-200:100:200]);yticks([-200:100:200]);
subplot(1,3,2);
plot(out_ang_exL4(a,3)*69-out_ang_exL4(a,1)*69,out_ang_inL4(a,3)*69-out_ang_inL4(a,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);
ref.Color='k';box off;xlim([-200 200]);ylim([-200 200]);hold on;title('L4');xticks([-200:100:200]);yticks([-200:100:200]);
subplot(1,3,3);
plot(out_ang_exL5(a,3)*69-out_ang_exL5(a,1)*69,out_ang_inL5(a,3)*69-out_ang_inL5(a,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',12);;
ref.Color='k';box off;xlim([-400 400]);ylim([-400 400]);hold on;title('L5');xticks([-400:200:400]);yticks([-400:200:400]);
%% Test
par1=abs(out_ang_exL23(a,3)*69-out_ang_exL23(a,1)*69);
par2=abs(out_ang_inL23(a,3)*69-out_ang_inL23(a,1)*69);
find(~isnan(par2))
P = polyfit(par1(find(~isnan(par2))),par2(find(~isnan(par2))),1);
    yfit = P(1)*par1(find(~isnan(par2)))+P(2);
 yfit2=1*par1(find(~isnan(par2)));
 
 [r t]=ttest2(yfit,yfit2)




%% Correlations Morphology and in vivo
com=[];com=[od_out_iviv(:,[1 2 3 4 5 7 8]) noise_corr' morph_parameters(:,1:21) sum_densap' sum_densba' max_densba' max_densap' MLba_diff' MLap_diff']
correlation_matrix(com,0);title('PMorphology od out');
%% 
m_res=[];
for i=1:length(nan_vector)
    if ~isnan(str(nan_vector(i)).morph(1))==1 & str(nan_vector(i)).resp==1 & ~isnan(str(nan_vector(i)).OSIpref)==1
        m_res(i)=1;
     
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=find(m_res==1);
%% 
%Ori and apical 
corr_plot(morph_parameters(morph_res_sub,3),od_out_iviv(morph_res_sub,4),pia_input,{'a','Apical width/height','Pial depth (µm)'});%xlim([-55 95]);ylim([0 180]);ylabel('Orientation preference','Color','k');xlabel('Angle centroid L2/3 IN','Color','b')
%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% Plot tuning curve for cell with low and high OSI where morphology is present
%pssible examples: 77,94 112,143
%Cell 143 %ipsi
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp15006\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex1=Fit(1).ipsi.FittedDataOri
ind_tr1=peaks(1).ipsi;
%% 
% load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp14373\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
% fit_ex1=Fit(1).contra.FittedDataOri
%% 
%Cell 122 contra
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp14757\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex2=Fit(2).contra.FittedDataOri;
ind_tr2=peaks(2).contra;
%% 
shift_ori(fit_ex1,fit_ex2,ind_tr1,ind_tr2,od_out_iviv(143,1),od_out_iviv(122,1),od_out_iviv(143,7),od_out_iviv(122,7));
%% 
id_m=143;
plot_morphologies(str,id_m,1,1);
%% 
%122
id_m=122;
plot_morphologies(str,id_m,1,1);
%% OSI and morph apical
corr_plot(morph_parameters(morph_res_sub,4),od_out_iviv(morph_res_sub,1),[],{'a','Apical width/height','Pial depth (µm)'});xlim([0 25]);ylim([0 1]);yticks([0:0.2:1]);ylabel('OSI','Color','k');xlabel('Number of branch points apical','Color','k');set(gca,'FontSize',12)
%% OSI and morph basal
corr_plot(morph_parameters(morph_res_sub,14),od_out_iviv(morph_res_sub,1),[],{'a','Apical width/height','Pial depth (µm)'});xlim([0 40]);ylim([0 1]);yticks([0:0.2:1]);ylabel('OSI','Color','k');xlabel('Number of branch points basal','Color','k');set(gca,'FontSize',12)
%% 
%% ODI and max branch order
morph_res_sub2=morph_res_sub;
find(od_out_iviv(morph_res_sub2,7)>50)
morph_res_sub2(10)=[]
corr_plot(morph_parameters(morph_res_sub2,2),od_out_iviv(morph_res_sub2,7),[],{'a','Apical width/height','Pial depth (µm)'});%xlim([0 25]);ylim([0 1]);yticks([0:0.2:1]);ylabel('OSI','Color','k');xlabel('Number of branch points apical','Color','k');set(gca,'FontSize',12)
%% Get all 147 morphtraces
for i=1:length(nan_vector)
    if ~isempty(str(nan_vector(i)).morph)==1 
    zz{:,i}=str(nan_vector(i)).morphtraces;
    else
        zz{:,i}=NaN;
    end
end
%% %% Morphology sholl analysis
close all;
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147);
%% Tuning width vs peak number of sholl crossings for apical
morph_res_sub2=morph_res_sub;
find(od_out_iviv(morph_res_sub2,7)>50)
morph_res_sub2(10)=[]
corr_plot(max_s(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'a','Apical width/height','Pial depth (µm)'});xlim([0 15]);;ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',12)
ylim([5 40]);yticks([0:10:40])
%% 
morph_res_sub2=morph_res_sub;
find(od_out_iviv(morph_res_sub2,7)>50)
morph_res_sub2(10)=[]
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'a','Apical width/height','Pial depth (µm)'});xlim([0 30]);;ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',12)
ylim([5 40]);yticks([0:10:40])

%% Sholl analysis for example cells 
temp=zz{143}
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(temp{1,1}, 20, '-s');
[sb, ddb, sdb, XP, YP, ZP, iD] = sholl_tree(temp{1,2}, 20, '-s');
temp=zz{122}
[s1, dd1, sd1, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,1}, 20, '-s');
[s1b, dd1b, sd1b, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,2}, 20, '-s');
figure;plot(dd,s,'-k');box off;
hold on;plot(ddb,sb,'--k')
hold on;plot(dd1,s1,'-b');
hold on;plot(dd1b,s1b,'--b');
legend('Apical','Basal');legend boxoff;
ylabel('Number of dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',12)
%% 
discretize_plot(pia_input,3,od_out_iviv(:,4),1);xlabel('Pial depth bins');ylabel('Orientation preference');xticks([1:1:3]);title('in vivo in vitro cells');ylim([0 180])
%% 
corr_plot(max_s(morph_cells_id)',L23fr(morph_cells_id,1),[],{'Number of branch points basal','Orientation selectivity index (gOSI)','Pial depth'});% xlim([0 22]);ylim([0 1])
%corr_plot(max_s(morph_cells_id)',sum(abs_inh_n(morph_cells_id,1:5),2),[],{'Number of branch points basal','Orientation selectivity index (gOSI)','Pial depth'});% xlim([0 22]);ylim([0 1])
%corr_plot(corr_exc_apical(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'Number of branch points basal','Orientation selectivity index (gOSI)','Pial depth'});% xlim([0 22]);ylim([0 1])
%% 
com=[morph_parameters(:,2) nanmax(df,[],2) dis_s'  morph_parameters(:,4)  max_s'  sum_densap' ...
    morph_parameters(:,12) nanmax(db,[],2) dis_s_ba'  morph_parameters(:,14) max_s_ba'  sum_densba' L23fr L4fr tot_input tot_inputL23 tot_inputL4]
   correlation_matrix(com,0)
   %% 
    corr_plot(nanmax(df(morph_cells_id,:),[],2),sum(frac_exv_n(morph_cells_id,1:5),2),[],{'Number of branch points basal','Orientation selectivity index (gOSI)',''}) 
   %% 
   corr_plot(nanmax(df(morph_cells_id,:),[],2),tot_input(morph_cells_id,1),[],{'Number of branch points basal','Orientation selectivity index (gOSI)',''}) 
   %% 
    corr_plot(morph_parameters(morph_cells_id,4),tot_input(morph_cells_id,1),[],{'Number of branch points apical','Mean EX input',''});ylabel('Mean EX input','Color','r')
    ylim([0 0.8]);yticks([0:0.2:0.8]);set(gca,'FontSize',12)
    %% 
corr_plot(morph_parameters(morph_cells_id,4),sum(frac_exv_n(morph_cells_id,1:10),2),[],{'Number of branch points basal','Orientation selectivity index (gOSI)',''}) 
    %% 

   a=find(od_out_iviv(:,1)<1)

com=[od_out_iviv(a,7) L23fr(a,:) L4fr(a,:) ex_spanhL23(a)' in_spanhL23(a)' ex_spanhL4(a)' in_spanhL4(a)' corr_exc_apical(a)']
correlation_matrix(com,0)
%% 
g1=[];
g2=[];
ia=morph_res_sub;  
%par=corr_inh_apical(ia)

g1=find(od_out_iviv(ia,1)<0.5)
% tmp=zeros(length(ia),1)
% tmp(g1)=g1;
g2=find(od_out_iviv(ia,1)>0.5)
%g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<360) 
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2])
xticklabels({'0-90°','90-180°'})
xlabel('Orientation');
%ylabel('|Angle Centroid L23|','Color','r');  set(gca,'FontSize',12);
%ylabel('|\Delta Centroid L23 X|','Color','r');  set(gca,'FontSize',12);
%ylabel('EX L4 fraction','Color','r');
%ylabel('IN L4 fraction','Color','b');
ylabel('EX L2/3 fraction','Color','r');
%ylabel('Pial depth (µm)')
set(gca,'FontSize',12);

%% 

%% Plot Morphologies of cells all or defined subgroups directly from structure
id_m=nan_vector(morph_cells_id)
plot_morphologies(str,id_m,10,10);


%% Plot morphologies of specific groups
plot_morphologies_iviv(zz,find(idx_input==1),10,10,m_flip_a)
%% Sort morpho by pial depth
[p spia]=sort(pia_input(morph_cells_id))
plot_morphologies_iviv(zz,morph_cells_id(spia),10,10,m_flip_a)
%% Plot only in vivo morpho cells
[ina inb]=intersect(id_m,find(iviv_cells==1));
plot_morphologies(str,ina,6,6);
ina=[];inb=[];
%% 
for i=1:156
    if ~isempty(str(i).morph)==1 & str(i).resp==1 & ~isempty(str(i).OSI)==1
        m_res(i)=1;
     
    else
       m_res(i)=NaN;
    end
end
morph_res=find(m_res==1);
%% 


%% 

%% Plot cells sorted based on OSI
plot_morphologies_iviv(zz,morph_res_sub,6,6)
%% 
[isa isb]=sort(od_out_iviv(morph_res_sub,1));
%% 
plot_morphologies_iviv(zz,morph_res_sub(isb),6,6,m_flip_a)
%% 
plot_morphologies_iviv(zz,[123 135],6,6,m_flip_a)


%% 
ylabel('Number of dendritic crossing');xlabel('Distance from Soma (µm)');set(gca,'FontSize',12);
%% Correlations Morphology and in vivo
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[morph_parameters(:,2) nanmax(df,[],2) dis_s'  morph_parameters(:,4)  max_s'  sum_densap' ...
    morph_parameters(:,12) nanmax(db,[],2) dis_s_ba'  morph_parameters(:,14) max_s_ba'  sum_densba'  od_out_iviv(:,[1 2 3 4 5 7 8]) pia_input]
G=correlation_matrix(com,0);
Gf=G(13:end,1:12)
Gf(6,1)=0;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 325]);imagesc(Gf);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:8]);
 
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Density','Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Density'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}','Pial depth'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% 



%% 

%% 
plot_morphologies_iviv(zz,[123 135],6,6,m_flip_a);
%% 
corr_plot(morph_parameters(id_ori,3),out_ang_inL23(id_ori,5),pia_input(id_ori),{'a','Apical width/height','Pial depth (µm)'})

%% Max sholl crossings apical
corr_plot(max_s(morph_res_sub)',od_out_iviv(morph_res_sub,1),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});ylabel('gOSI','Color','k');xlabel('Peak number of sholl crossings apical','Color','k');set(gca,'FontSize',12);yticks(0:0.25:1)
%% Max sholl crossings basal
corr_plot(max_s_ba(morph_res_sub)',od_out_iviv(morph_res_sub,1),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});ylabel('gOSI','Color','k');xlabel('Peak number of sholl crossings basal','Color','k');set(gca,'FontSize',12);yticks(0:0.25:1)
%% Tuning width sholl corssings apical
corr_plot(max_s(morph_res_sub)',od_out_iviv(morph_res_sub,7),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});ylabel('Tuning width (°)','Color','k');xlabel('Peak number of sholl crossings apical','Color','k');set(gca,'FontSize',12);
%% Tuning width sholl corssings
corr_plot(max_s_ba(morph_res_sub)',od_out_iviv(morph_res_sub,7),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});ylabel('Tuning width (°)','Color','k');xlabel('Peak number of sholl crossings basal','Color','k');set(gca,'FontSize',12);
%% 
corr_plot(max_s(morph_res_sub)',od_out_iviv(morph_res_sub,7),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});
%% Plot individual cell with circles
a=122;
ml=13;
  tmp=zz{1,a}
  figure;set(gcf,'color','w')
  m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
         m.EdgeColor = [0 0 1]
        m1=plot_tree(tmp{1,2},[0 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
        m1.EdgeColor = [0.5 0.5 0.5]
m2=plot_tree(tmp{1,3},[0 0 1],[0 tmp{1,5} 0],[],1,'-b');hold on;
        m2.EdgeColor = [0.8 0.8 0.8]
set(gca,'Ydir','reverse');
 axis off
%% 

for i=1:ml
hold on;
  viscircles([0 pia_input(a)],20*i,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
 axis off
 %% Histogram 
 figure;set(gcf,'color','w')
 h4 = histogram(morph_parameters(:,4),10);box off;h4.EdgeColor = 'k';h4.FaceColor = 'k';xlabel('');ylabel('Counts')
hold on;h4 = histogram(morph_parameters(:,14),10);box off;h4.EdgeColor = 'k';h4.FaceColor = 'm';xlabel('Number of branch points');ylabel('Counts')
 legend('Apical','Basal');set(gca,'FontSize',12);legend boxoff;set(gca,'FontSize',12);


%% 


%% 


com=[];com=[L23fr(iaa,1) L4fr(iaa,1)  L23fr(iaa,2) L4fr(iaa,2) out_ang_exL23(iaa,5)  out_ang_exL4(iaa,5) out_ang_inL23(iaa,5)  out_ang_inL4(iaa,5) out_ang_exL23(iaa,3)  out_ang_exL4(iaa,3) out_ang_inL23(iaa,3)  out_ang_inL4(iaa,3) od_out_iviv(iaa,[1 2 3 4 5 6 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:20]);yticks([1:1:20]);
%% 
figure;histogram(od_out_iviv(iaa,4))

%% Important correlations
corr_plot(out_ang_inL23(:,5),od_out_iviv(:,4),pia_input,{'a','Apical width/height','Pial depth (µm)'});xlim([-55 95]);ylim([0 180]);ylabel('Orientation preference','Color','k');xlabel('Angle centroid L2/3 IN','Color','b')
%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% 
corr_plot(od_out_iviv(id_ori,4),out_ang_inL23(id_ori,5),idx_input(id_ori),{'a','Apical width/height','Pial depth (µm)'});ylim([-60 100]);xlim([0 180]);xlabel('Orientation preference','Color','k');ylabel('Angle centroid L2/3 IN','Color','b')
set(gca,'FontSize',12);%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% 
corr_plot(od_out_iviv(id_ori,4),out_ang_exL4(id_ori,5),pia_input(id_ori),{'a','Apical width/height','Pial depth (µm)'});ylim([45 95]);xlim([0 180]);xlabel('Orientation preference','Color','k');ylabel('Angle centroid L4 EX','Color','r')
set(gca,'FontSize',12);%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% lower and upper part of pia

corr_plot(od_out_iviv(find(pia_input>100),4),L4fr(find(pia_input>100),1),pia_input(find(pia_input>100)),{'a','Apical width/height','Pial depth (µm)'});
%% 
corr_plot(od_out_iviv(find(pia_input<240),4),out_ang_exL4(find(pia_input<240),5),pia_input(find(pia_input<240)),{'a','Apical width/height','Pial depth (µm)'});ylim([45 95]);xlim([0 180]);xlabel('Orientation preference','Color','k');ylabel('Angle centroid L2/3 IN','Color','b')
%% 
corr_plot(od_out_iviv(find(pia_input>220),4),pia_input(find(pia_input>220)),[],{'a','Apical width/height','Pial depth (µm)'});
%% 
a=300
corr_plot(od_out_iviv(find(pia_input<a),4),pia_input(find(pia_input<a)),gd(find(pia_input<a)),{'a','Apical width/height','Pial depth (µm)'});

%% 
corr_plot(od_out_iviv(id_ori,4),out_ang_inL23(id_ori,3),gd(id_ori),{'a','Apical width/height','Pial depth (µm)'});
set(gca,'FontSize',12);%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% %%%%%%%%CLUSTERING

%% %% Clustering HCA VERTICAL
clu_num =4;
%pcs =[];
pcs     =[1 2 3];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input, clustering_input, leafOrder] = hca([score_ex(:,pcs) score_in(:,pcs)],0,'ward',clu_num,pia_input,2,0.6);%call function for clustering
%[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%% 
clu_num =5;
%pcs =[];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input, clustering_input, leafOrder] = hca([L23fr L4fr L5fr],0,'ward',clu_num,pia_input,0,0.6);%call function for clustering
%[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%% 



%% Plot average maps of clusters
plot_avg_maps(str,nan_vector,ex_map,in_map,pia_input,10,1,idx_input);

%% kmeans clustering
% idx_input=[];
% idx_input = kmeans([score_ex(:,pcs) score_in(:,pcs)],4)
%% Barplot difference of clusters
%Pial depth
[statsout] = barplot_sw(pia_input,idx_input,{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45)
%% Orientation preference
a=find(od_out_iviv(:,1)>0.25)  
[statsout] = barplot_sw(od_out_iviv(a,4),idx_input(a),{'Clusters','Orientation preference'});xtickangle(45);set(gca,'Ydir','reverse');set(gca,'FontSize',12)
%% Orientation selectivity
[statsout] = barplot_sw(od_out_iviv(:,6),idx_input,{'Clusters','Ca_{peak}'})
%% 
[statsout] = barplot_sw(noise_corr',idx_input,{'Clusters','Ca_{peak}'})
%% Angle centroid
[statsout] = barplot_sw(abs(out_ang_inL23(:,3)-out_ang_inL23(:,1)),idx_input,{'Clusters','Angle Centroid L23in'})
%% L4 ex fraction
[statsout] = barplot_sw(L5fr(:,1),idx_input,{'Clusters','L4 ex fraction'});set(gca,'FontSize',12)
%% 
%% 

[statsout] = barplot_sw(in_spanhL23',idx_input,{'Clusters','L4 ex fraction'});set(gca,'FontSize',12)
%% L4 total input sum
[statsout] = barplot_sw(abs(tot_inputL4(:,1)),idx_input,{'Clusters','L4 total input'})
%% 
[statsout] = barplot_sw(sum(abs_exv_n(:,6:7),2),idx_input,{'Clusters','L4 total input'})

%% L23 in fraction
[statsout] = barplot_sw(L4fr(:,1),idx_input,{'Clusters','L23 in fraction'})
%% L23 in fraction
[statsout] = barplot_sw(tot_inputL23(:,1),idx_input,{'Clusters','L23 total input'})
%% 
[statsout] = barplot_sw(tot_input,idx_input,{'Clusters','Angle IN L23'})
%% 
[statsout] = barplot_sw(ex_spanh',idx_input,{'Clusters','Angle IN L23'})
%% 
[statsout] = barplot_sw(out_ang_exL4(id_ori,5),idx_input(id_ori),{'Clusters','L23 total input'})
%% 
[statsout] = barplot_sw(pia_input(id_ori),idx_input(id_ori),{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45);ylim([0 400])
%% 
[statsout] = barplot_sw(L4fr(id_ori,1),idx_input(id_ori),{'Clusters','L4 ex fraction'});
%% 
[fa]=discretize(pia_input(id_ori),3)
[statsout] = barplot_sw(od_out_iviv(id_ori,4),fa,{'Pial depth bins','Orientation preference'});xtickangle(45);set(gca,'Ydir','reverse');set(gca,'FontSize',12)
%% 
[fa]=discretize(pia_all,2)
[statsout] = barplot_sw(od_out(:,3),fa',{'Pial depth bins','Orientation preference'});xtickangle(45);set(gca,'Ydir','reverse');set(gca,'FontSize',12)
%% 
%% 
[statsout] = barplot_sw(od_out_iviv(a,4),idx_input(a),{'Pial depth bins','Orientation preference'})
%% 

[statsout] = barplot_sw(ex_spanhL5',idx_input,{'Pial depth bins','Orientation preference'})
%% 
[statsout] = barplot_sw(L23fr(id_ori,1),idx_input(id_ori),{'Clusters','L23 in fraction'})
%%  morphology
[statsout] = barplot_sw(morph_parameters(:,21),idx_input,{'Clusters','L23 total input'})
%% 
[statsout] = barplot_sw(ap_corrin',idx_input,{'Clusters','Correlation morpho input'})
%% 
[statsout] = barplot_sw(morph_parameters(id_ori,3),idx_input(id_ori),{'Clusters','Correlation morpho input'})
%% Plot clusters with values underneath it as heatmaps
dendroplot_SW(clustering_input,leafOrder,0.8,[score_ex(:,1:3) score_in(:,1:3)],{'PC1ex','PC2ex','PC3ex','PC1in','PC2in','PC3in','Slice Ori','setups'},pia_input)
%% Looking at input fractions
dendroplot_SW(clustering_input,leafOrder,0.8,[L23fr(:,1) L4fr(:,1)  L5fr(:,1) L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3ex','L4ex','L5ex','L2/3in','L4in','L5in'},pia_input)
%% Looking at input fractions ex
dendroplot_SW(clustering_input,leafOrder,11,[L23fr(:,1) L4fr(:,1)  L5fr(:,1)],{'L2/3ex','L4ex','L5ex'})
%% Looking at input fractions in
dendroplot_SW(clustering_input,leafOrder,11,[L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3in','L4in','L5in', 'slice_ori','setups'})
%% Plot clusters with values underneath it as heatmaps: SETUPS
dendroplot_SW(clustering_input,leafOrder,0.8,[od_out_iviv(:,4)],{'ORI'})
%% 
%% Plot clusters with values underneath it as heatmaps: SETUPS
dendroplot_SW(clustering_input,leafOrder,0.8,[setups],{'Oripref'})

%% Plot clusters with values underneath it as heatmaps: SLICE ORIENTATION

dendroplot_SW(clustering_input,leafOrder,11,[slice_ori' score_in(:,3)],{'Slice Ori','PC3in'})
%% 
com=[];com=[ap_correx' ap_corrin' sftf_out_iviv];
correlation_matrix(com,0);title('correlation input morpho od_out iviv');
%% 
corr_plot(ap_correx',sftf_out_iviv(:,4),pia_input,{'L4ex input','Apical width/height','XSA'}); 
%% Histograms for pial depth bins for all in vivo and in vivo invitro
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 200]);
subplot(1,2,1);
[gd] =  discretize(pia_all,3);
binRange = 0:45:180;
hcx = histcounts(od_out(find(gd==1),3),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out(find(gd==2),3),[binRange Inf],'Normalization','probability');
hcz = histcounts(od_out(find(gd==3),3),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy;hcz]');box off;xlabel('Orientation (deg)','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
b1(3).FaceColor='m';
ylim([0 0.2]);yticks([0:0.1:0.2]);ylabel('Fraction of cells');set(gca,'FontSize',10);hold on;title('all cells')
subplot(1,2,2);
gd=[];
[gd] =  discretize(pia_input,3);
hcx = histcounts(od_out_iviv(find(gd==1),4),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(gd==2),4),[binRange Inf],'Normalization','probability');
hcz = histcounts(od_out_iviv(find(gd==3),4),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy;hcz]');box off;xlabel('Orientation (deg)','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
b1(3).FaceColor='m';
set(gca,'FontSize',10)
%ylim([0 0.2]);yticks([0:0.1:0.2]);hold on;title('iviv cells');legend('upper L2/3','mid L2/3','lower L2/3');legend boxoff
%% 


  

%% Morphology and in vivo

 %% 
 
 