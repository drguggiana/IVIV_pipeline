%Load structure 
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
    bin_exc = str(nan_vector(cells)).excMap(3:16,:);
    bin_excm= [zeros(2,16); bin_exc];
    bin_inh = str(nan_vector(cells)).inhMap;
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
    bin_exc = str(nan_vector(cells)).excMap(3:16,:);
    bin_excm= [zeros(2,16); bin_exc]<0;
    bin_inh = str(nan_vector(cells)).inhMap>0;
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
    if ~isempty(str(nan_vector(i)).morph)==1;
        morph_cells(i)=1;
        m_flip_a(i)=str(nan_vector(i)).morph_flip_again;
        morph_parameters(i,:)=str(nan_vector(i)).morph;
    else
        morph_cells(i)=0;
        m_flip_a(i)=NaN;
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
L5fr(find(L5fr(:,1)==0),1)=NaN ;
L5fr(find(L5fr(:,2)==0),2)=NaN ;
diffL23fr=L23fr(:,1)-L23fr(:,2);
diffL4fr=L4fr(:,1)-L4fr(:,2);
diffL5fr=L5fr(:,1)-L5fr(:,2);
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
tmp1=ex_map_raw(:,:,i);
tmp2=in_map_raw(:,:,i);
% tmp1=ex_map(:,:,i);
% tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
%tot_input(i,:)=[sum(tmp1(:))/length(tmp1);sum(tmp2(:))/length(tmp2)];
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map_raw(3:5,:,i);
tmp2=in_map_raw(3:5,:,i);
%  tmp1=ex_map(3:5,:,i);
%  tmp2=in_map(3:5,:,i);
tot_inputL23(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map_raw(6:7,:,i);
tmp2=in_map_raw(6:7,:,i);
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
hold on;plot(nanmedian(diffL23fr),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');
subplot(3,2,3);
h4 = histogram(diffL4fr,8);h4.BinWidth =  0.025;h4.EdgeColor = 'k';h4.FaceColor = 'g'
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylabel('Cell counts');hold on;set(gca,'FontSize',10);box off;
text(-0.2,40,'L4');
hold on;plot(nanmedian(diffL4fr),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
xlim([-0.3 0.3]);
subplot(3,2,5);
h4 = histogram(nonzeros(diffL5fr),8);h4.BinWidth =  0.025;h4.FaceColor = [0.5 0.5 0.5];ylabel('Cell counts');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');
xlim([-0.3 0.3]);xlabel('\Delta EX-IN vertical fraction');hold on;set(gca,'FontSize',10);
text(-0.2,40,'L5','Color','k');
hold on;plot(nanmedian(diffL5fr),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
%Horizontal span
subplot(3,2,2);
h4 = histogram((ex_spanhL23-in_spanhL23)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'm';
xlim([-10*xf 10*xf]);hold on; title('Horizontal');box off;
hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,50,'IN','Color','b');text(7*xf,50,'EX','Color','r');
text(-7*xf,40,'L2/3','Color','k');
hold on;plot(nanmedian((ex_spanhL23-in_spanhL23)*xf),50,'v','MarkerFaceColor','m','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,4);
hold on;h4 = histogram((ex_spanhL4-in_spanhL4)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = 'g';
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7*xf,40,'L4','Color','k');
hold on;plot(nanmedian((ex_spanhL4-in_spanhL4)*xf),50,'v','MarkerFaceColor','g','MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
subplot(3,2,6);
hold on;h4 = histogram((ex_spanhL5-in_spanhL5)*xf,8);h4.BinWidth = 1*xf;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];;
xlim([-10*xf 10*xf]);hold on;line([0 0], [1 50],'Color','k','LineStyle','--');set(gca,'FontSize',10);
xlabel('\Delta EX-IN horizontal span (µm)');hold on;set(gca,'FontSize',10);
text(-7*xf,40,'L5','Color','k');
hold on;plot(nanmedian((ex_spanhL5-in_spanhL5)*xf),50,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');ylim([0 50]);yticks([0:25:50])
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
%% Calculate centroid of map and angle towards it 
for i=1:length(nan_vector)
somax(i)=str(nan_vector(i)).somaCenter(1);
somay(i)=552-str(nan_vector(i)).somaCenter(2);
end
%fake map
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
ang1=out_ang_exL4;
ang2=out_ang_inL4;
idxtp=1:147;
%call function quiver
quiver_centroid(ang1,ang2,idxtp,ex_map,in_map,119);
%% Display soma posiiton and how it was aligned for PCA
idxtp=1:147;
ang1=out_ang_exL23;
ang2=out_ang_inL23;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 900, 200]);
subplot(1,3,1);
for i=1:length(idxtp)
 hold on;plot(ang1(i,1),ang1(i,2),'k^','MarkerFaceColor','w','MarkerSize',5);   
end
hold on;xlim([4 12]);ylim([1 8]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8 8],'Color','k','LineStyle','--');hold on;line([8 8], [1 16],'Color','k','LineStyle','--');axis off;box off;set(gca,'Ydir','reverse');
subplot(1,3,2);
for i=1:length(idxtp)
 hold on;plot(8,ang1(i,2),'k^','MarkerFaceColor','w','MarkerSize',5);   
end
hold on;xlim([4 12]);ylim([1 8]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8 8],'Color','k','LineStyle','--');hold on;line([8 8], [1 16],'Color','k','LineStyle','--');axis off;box off;set(gca,'Ydir','reverse');
subplot(1,3,3);
for i=1:length(idxtp)
 hold on;plot(ang1(i,1),4,'k^','MarkerFaceColor','w','MarkerSize',5);   
end
hold on;xlim([4 12]);ylim([1 8]);hold on;line([1 16], [2 2],'Color','k','LineStyle','--');hold on;line([1 16], [6 6],'Color','k','LineStyle','--');
hold on;line([1 16], [8 8],'Color','k','LineStyle','--');hold on;line([8 8], [1 16],'Color','k','LineStyle','--');axis off;box off;set(gca,'Ydir','reverse');
%% Distributions of actual centroid all input cells for EX and IN
a=1:147;
centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,0,[],{});
%% % Relation betwen ex and in centroid
ang1=out_ang_exL4;
ang2=out_ang_inL4;
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
%% Desired correlation between angle and inhbition PC2in
corr_plot(out_ang_inL23(:,5),score_in(:,2),pia_input,{'Angle centroid L2/3 IN','PC2_{in}','Pial depth (µm)'});ylabel('PC2_{in}','Color','b');xlabel('Angle centroid L2/3 IN','Color','b');
%% Desired correlation between horizontal centroid and PC3in
corr_plot((out_ang_inL23(:,3)-out_ang_inL23(:,1))*69,score_in(:,3),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});xlim([-150 150]);ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% 
corr_plot(sum(abs_exv_n(:,1:5),2),score_ex(:,1),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});;ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% 
corr_plot(sum(abs_inv_n(:,1:3),2),score_in(:,1),pia_input,{'Horzontal centroid L2/3 IN','PC3_{in}','Pial depth'});;ylabel('PC3_{in}','Color','b');xlabel('Horzontal centroid L2/3 IN','Color','b');
%% FAKE MAPS CENTROID
corr_plot(out_ang_fake(:,5),score_in(:,2),pia_input,{'Angle centroid L23 fake','PCin2','Pial depth'});set(gca,'Ydir','reverse');
%% 

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
%% %% Morphology parameters vs real parameters vertical
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) tot_inputL23(:,1) tot_inputL4(:,1) tot_inputL23(:,2) tot_inputL4(:,2) out_ang_inL23(:,5) out_ang_inL23(:,3) morph_parameters(:,1:21) sum_densap' sum_densba']; 
correlation_matrix(com,0);title('Morphology vs Vertical input');
%% Total input apical dendrite branch points
corr_plot(morph_parameters(:,4),abs(tot_input(:,1)),pia_input,{'Branch points','Total input EX','Pial depth (µm)'})
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
%% Correlations PCs input and OD out
ia=find(od_out_iviv(:,1)>0.3)  
%b=find(od_out_iviv(:,1)>0.3)  
%[ia ib]=intersect(a,b)
com=[];com=[score_ex(ia,1:3) score_in(ia,1:3) od_out_iviv(ia,:)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% 
%% Correlations PCs input and OD out
a=find(od_out_iviv(:,2)>0.3)  
b=find(od_out_iviv(:,1)>0.3)  
[ia ib]=intersect(a,b)
com=[];com=[score_com(ia,1:10) od_out_iviv(ia,:)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% 
a=find(od_out_iviv(:,1)>0.3)  
 %b=find(L4fr(:,1)>0)  
 %[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(a,4),sum(abs_exv_m(a,6),2),pia_input(a),{'Direction','PC3ex','pia'})
%% 
a=find(od_out_iviv(:,1)>0.3)  
b=find(sum(abs_inv(:,1:16),2)<9000)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(ia,1),sum(abs_inv(ia,1:16),2),pia_input(ia),{'Direction','PC3ex','pia'})
%% 
corr_plot(od_out_iviv(ia,5),score_ex(ia,3),pia_input(ia),{'Direction','PC3ex','pia'})
%% 
corr_plot(od_out_iviv(ia,5),frh_medial_s(ia,1),pia_input(ia),{'Direction','Medial Fraction EX','pia'})

%% Correlations PCs input and sftft out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) sftf_out_iviv(:,2:end)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations Input vertical vs ODout
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) tot_input(:,1) tot_input(:,2) out_ang_inL23(:,3) od_out_iviv(:,1:8)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% 
a=find(od_out_iviv(:,1)>0.3)  
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) tot_input(a,1) tot_input(a,2) out_ang_inL23(a,3)  sum(abs_exv_n(a,1:5),2) ex_spanhL23(a)' ex_spanhL4(a)' frh_lateral_s(a,1) frh_lateral_s(a,2) od_out_iviv(a,1:8)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Score in with Orientation
a=find(od_out_iviv(:,1)>0.3)  
b=find(pia_input<200)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(ia,4),score_in(ia,2),pia_input(ia),{'11','2','pia'})
corr_plot(od_out_iviv(ia,4),out_ang_inL23(ia,5),pia_input(ia),{'Orientation','Angle centre of mass L23 in','Pial depth'})
;set(gca,'FontSize',12);
%% OSI against total input
a=find(od_out_iviv(:,1)>0.3)  
b=find(tot_input(:,2)<200)
[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(ia,1),tot_input(ia,2),[],{'OSI','IN INPUT','pia'});

%% 
a=find(od_out_iviv(:,1)>0.3)  
corr_plot(od_out_iviv(a,5),out_ang_exL23(a,3)-out_ang_inL23(a,3),pia_input(a),{'OSI','IN INPUT','pia'});
%% Correlations Input vertical vs ODout
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) tot_input(:,1) od_out_iviv]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% 
com=[];com=[ex_spanhL23' in_spanhL23' ex_spanhL4' in_spanhL4'  abs(ex_spanhL23'-in_spanhL23')*69 od_out_iviv]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');

%% Correlations Input vertical vs SFTF pref
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) sftf_out_pref_iviv(:,1:3)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% 



%% Correlations Input vertical vs ODout
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) tot_input(a) tot_inputL23(a) tot_inputL4(a) frh_lateral_s(a,1) frh_lateral_s(a,2) od_out_iviv(a,1:8)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');

%% Correlations Input vertical vs ODout
com=[];com=[L23h L4h L5h od_out_iviv(:,:)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');

%% 
corr_plot(sftf_out_iviv(:,2),abs(out_ang_exL23(:,3)-out_ang_inL23(:,3)),[],{'SF','Centroid L23 EX-IN','XSA'});
 
%% 
a=find(od_out_iviv(:,2)>0.25)
corr_plot(od_out_iviv(a,5),in_spanhL23(a)',[],{'SF','Centroid L23 EX-IN','XSA'});

%% SF vs abs difference between L23 EX/IN span
%par=abs(ex_spanhL23-in_spanhL23)*69
par=max_densap
%par=score_ex(:,1)
% g1=find(sftf_out_iviv(:,2)==0.02)
% g2=find(sftf_out_iviv(:,2)==0.08)
 g1=find(od_out_iviv(:,1)<0.5)
 g2=find(od_out_iviv(:,1)>0.5)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'0.02','0.08'})
;xlabel('SF');ylabel('\Delta EX-IN horizontal span (µm)');
%ylim([0 200]);
%% Correlations Input horizontal vs ODout
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' od_out_iviv]; 
correlation_matrix(com,0);title('Input horizontal ODout iviv');
%% Correlations PCs input and spon out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) spon_out_iviv]; 
correlation_matrix(com,1);title('PC inputs spon_out iviv');
%% Correlations Vertical input and spon out
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov spon_out_iviv]; 
correlation_matrix(com,1);title('Input vertical spon_out iviv');
%% Correlations Horizontal input and spon out
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)'  spon_out_iviv]; 
%% Correlations PCs input and sftf_out_iviv
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov sftf_out_iviv(:,[1:6])]; 
correlation_matrix(com,0);title('PC inputs sftf_out iviv');
%% Correlations PCs input and sftf_out_iviv
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' sftf_out_iviv]; 
correlation_matrix(com,1);title('PC inputs sftf_out iviv');
%% Check correlation between ORI Pref with OSI cutoff
a=find(od_out_iviv(:,1)>0.3)
corr_plot(od_out_iviv(a,4),L4fr(a,1),[],{'ORI','L4ex fraction'});set(gca,'FontSize',12);
%% %% Intersectional approach ORI
a=find(L4fr(:,1)<0.3)  
b=find(od_out_iviv(:,1)>0.3)
c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
[ia ib]=intersect(a,b)
[iaa ic]=intersect(ia,c);
corr_plot(od_out_iviv(iaa,4),out_ang_inL23(iaa,5),[],{'ORI','L4ex fraction'});set(gca,'FontSize',12);
%% Correlation matrix with intersectional approach 
com=[];com=[score_ex(iaa,1:3) score_in(iaa,1:3) od_out_iviv(iaa,[1 2 3 4 5 6 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:20]);yticks([1:1:20]);
com=[];com=[L23fr(iaa,1) L4fr(iaa,1)  L23fr(iaa,2) L4fr(iaa,2) out_ang_exL23(iaa,5)  out_ang_exL4(iaa,5) out_ang_inL23(iaa,5)  out_ang_inL4(iaa,5) out_ang_exL23(iaa,3)  out_ang_exL4(iaa,3) out_ang_inL23(iaa,3)  out_ang_inL4(iaa,3) od_out_iviv(iaa,[1 2 3 4 5 6 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:20]);yticks([1:1:20]);

%%  Check correlation between DIR Pref with DSI cutoff
a=find(od_out_iviv(:,2)>0.3)
corr_plot(od_out_iviv(a,5),(out_ang_inL23(a,3)-out_ang_inL23(a,1))*69,slice_ori(a),{'DIR','Horizontal centre of mass','Slice Ori'});set(gca,'FontSize',12);
%% 
a=find(slice_ori>0)
corr_plot(od_out_iviv(a,4),L4fr(a,1),pia_input(a),{'DIR','Horizontal centre of mass','Slice Ori'});set(gca,'FontSize',12);
%% 

a=find(slice_ori>=0)  
b=find(od_out_iviv(:,1)>0.25)
[ia ib]=intersect(a,b)
corr_plot(od_out_iviv(ia,4),slice_ori(ia),pia_input(ia),{'DIR','Horizontal centre of mass','Slice Ori'});set(gca,'FontSize',12);
%% Check correlation between DIR Pref with DSI cutoff
a=find(od_out_iviv(:,2)>0.3)
corr_plot(od_out_iviv(a,5),frh_lateral_s(a),[],{'DIR','Fraction Lateral'});set(gca,'FontSize',12);
%% %% Intersectional approach ORI
a=find(od_out_iviv(:,2)>0.3)  
b=find((out_ang_inL23(:,3)-out_ang_inL23(:,1))*69>20)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
[ia ib]=intersect(a,b)
[iaa ic]=intersect(ia,c);
corr_plot(od_out_iviv(ia,5),(out_ang_inL23(ia,3)-out_ang_inL23(ia,1))*69,[],{'DIR','Horizontal centre of mass'});set(gca,'FontSize',12);
%% DIstributions of actual centroid
%a=find(od_out_iviv(:,2)>0.3);
%a=find(od_out_iviv(:,2)>0.3)
%a=1:147
a=find(od_out_iviv(:,1)>0.3)
%fe=od_out_iviv(a,4)
fe=slice_ori(a)
centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,1,fe,{'DIR'})
%% 

corr_plot(od_out_iviv(a,4),abs((out_ang_inL23(a,3)-out_ang_inL23(a,1)))*69,slice_ori(a),{'DIR','Horizontal centre of mass','slice_ori'});set(gca,'FontSize',12);

%% 
corr_plot(od_out_iviv(a,4),out_ang_inL23(a,5),pia_input(a),{'DIR','Horizontal centre of mass','slice_ori'});set(gca,'FontSize',12);
%% Check if orthognal orientations have different input arrangement
%Centroid in x and y
a=find(od_out_iviv(:,1)>0.25)
par=abs((out_ang_exL4(a,3)))*69;
%par=out_ang_exL4(a,7)
%par=abs(out_ang_inL23(a,8))*69;
%par=score_ex(:,2);
%par=pia_input(a)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) 
g2=find(od_out_iviv(a,4)>110 & od_out_iviv(a,4)<160) 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'~45°','~135'})
xlabel('Orientation');ylabel('L23 ex Centre of mass X (µm)');  set(gca,'FontSize',12);
%% 
%% Check if orthognal orientations have different input arrangement
%Centroid in x and y
a=find(od_out_iviv(:,1)>0.3)
%par=abs(out_ang_exL23(a,3)-out_ang_exL23(a,1))*69;
%par=abs(out_ang_exL23(a,5));
par=frh_lateral_s(a,1)
%par=out_ang_exL4(a,7)
%par=abs(out_ang_inL23(a,8))*69;
%par=score_ex(:,2);
%par=pia_input(a)
g1=find(od_out_iviv(a,4)>=0 & od_out_iviv(a,4)<30) 
g2=find(od_out_iviv(a,4)>80 & od_out_iviv(a,4)<110) 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'~0°','~90'})
xlabel('Orientation');ylabel('L23 ex COM X');  set(gca,'FontSize',12);

%% %% 
%Centroid in x and y
a=find(od_out_iviv(:,1)>0.25)
%par=abs((out_ang_exL23(a,4)-out_ang_exL23(a,2)))*69;
%par=L4h(a,3);
%par=abs(out_ang_inL23(a,8))*69;
par=score_in(a,2);
%par=pia_input(a)
g1=find(od_out_iviv(a,4)>0 & od_out_iviv(a,4)<90) 
g2=find(od_out_iviv(a,4)>90 & od_out_iviv(a,4)<180) 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'~45°','~135'})
xlabel('Orientation');ylabel('L23 ex Centre of mass X (µm)');  set(gca,'FontSize',12);
%% 
ia=find(od_out_iviv(:,1)>0.25)
%par=abs((out_ang_exL23(a,4)-out_ang_exL23(a,2)))*69;
%par=L4h(a,3);
%par=abs(out_ang_inL23(a,8))*69;
%par=abs((out_ang_exL4(a,4)-out_ang_exL4(a,2)))*69;
%par=abs(out_ang_exL4(a,3))
%par=out_ang_exL23(a,3)
%par=pia_input(a)
%par=L4fr(a,1)

%par=tot_input(a,2)
%b=find(pia_input>180 & pia_input<280)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
%[ia ib]=intersect(a,b)

%par=abs((out_ang_exL4(ia,3)-out_ang_exL4(ia,1)))*69;
%par=L4fr(ia,1)
par=sum(abs_exv_n(ia,6:7),2)
%par=out_ang_exL23(ia,5)
% g1=find(od_out_iviv(ia,4)>45 & od_out_iviv(ia,4)<90) 
% g2=find(od_out_iviv(ia,4)>135 & od_out_iviv(ia,4)<180) 
 g1=find(od_out_iviv(ia,4)>0 & od_out_iviv(ia,4)<90) 
 g2=find(od_out_iviv(ia,4)>90 & od_out_iviv(ia,4)<180) 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'~45°','~135'})
xlabel('Orientation');ylabel('L23 ex Centre of mass X (µm)');  set(gca,'FontSize',12);

%% EX L4 fraction to orientation
corr_plot(od_out_iviv(ia,4),sum(abs_inv_n(ia,1:5),2),pia_input(ia),{'Orientation','L4 EX fraction','ODI','Pia input'});set(gca,'FontSize',12);ylabel('L4 EX fraction','Color','r');xlabel('Orientation');xlim([0 180])
%% 
corr_plot(od_out_iviv(ia,4),abs(ex_spanhL4(ia)'-in_spanhL4(ia)'),pia_input(ia),{'Orientation','L4 EX fraction','ODI','Pia input'});set(gca,'FontSize',12);ylabel('L4 EX fraction','Color','r');xlabel('Orientation');xlim([0 180])

%% EX L4 total input to orientation
corr_plot(od_out_iviv(ia,4),abs(tot_inputL4(ia,1)),[],{'Orientation','L4 EX fraction','ODI'});set(gca,'FontSize',12);ylabel('L4 EX total input','Color','r');xlabel('Orientation');xlim([0 180]) 
%% EX L4 total input to orientation
corr_plot(od_out_iviv(ia,4),pia_input(ia),[],{'Orientation','L4 EX fraction','ODI'});set(gca,'FontSize',12);ylabel('L4 EX total input','Color','r');xlabel('Orientation');xlim([0 180])
%% 
ia=find(od_out_iviv(:,1)>0.25)  
%b=find(pia_input<210)
%c=find(od_out_iviv(:,4)<135)
%c=find(pia_input<300)
%[ia ib]=intersect(a,b)

corr_plot(spon_out_iviv(ia,2),ex_spanhL4(ia)',pia_input(ia),{'Orientation','L4 EX fraction','ODI','Pia input'});set(gca,'FontSize',12);ylabel('L4 EX fraction','Color','r');xlabel('Orientation');
%% Direction differences
a=find(od_out_iviv(:,2)>0.3)
%par=abs((out_ang_inL23(a,3)-out_ang_inL23(a,1)))*69;
%par=abs(out_ang_exL4(a,5)-out_ang_inL4(a,5));
%par=score_ex(a,3);
%par=abs(out_ang_inL23(a,5))
%par=score_ex(:,2);
%par=pia_input(a)
par=frh_medial_s(a,2)
 g1=find(od_out_iviv(a,5)>25 & od_out_iviv(a,5)<80) 
 g2=find(od_out_iviv(a,5)>250 & od_out_iviv(a,5)<360) 
 [statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
 xticklabels({'~45°','~315'})
xlabel('Direction');ylabel('Total in fraction medial');  set(gca,'FontSize',12);
%xlabel('Direction');ylabel('Pial depth (µm)');  set(gca,'FontSize',12);
%% 

par=abs(out_ang_inL23(a,4))*69;
g1=find(od_out_iviv(a,5)>45 & od_out_iviv(a,5)<225) 
tmp=zeros(length(a),1)
tmp(g1)=g1
g2=find(tmp==0)
%g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<360) 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'45°','~135'})
xlabel('Direction');ylabel('L23 ex Centre of mass X (µm)');  set(gca,'FontSize',12);


%% 

a=find(od_out_iviv(:,1)>0.25)
b=find(od_out_iviv(:,4)>20 & od_out_iviv(:,4)<70) 
c=find(od_out_iviv(:,4)>110 & od_out_iviv(:,4)<160) 
[ia ib]=intersect(a,b);
[ic id]=intersect(a,c);


corr_plot(od_out_iviv([ia' ic'],4),out_ang_inL23([ia' ic'],5),pia_input([ia' ic']),{'DIR','Horizontal centre of mass','slice_ori'});set(gca,'FontSize',12);
   %% 

%par=frh_lateral_s(:,1)
%par=abs((out_ang_inL4(a,3)-out_ang_inL4(a,1)))*69
%par=slice_ori(a)
%par=abs((out_ang_exL4(a,4)-out_ang_exL4(a,2)))*69;
%par=pia_input(a)
%par=score_in(:,3)
%par=in_spanhL23(a)*69
%   g1=find(od_out_iviv(a,5)>45 & od_out_iviv(a,5)<135) 
%   g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<300) 
%  g1=find(od_out_iviv(a,5)>45 & od_out_iviv(a,5)<225) 
%  g2=find(od_out_iviv(a,5)>225 & od_out_iviv(a,5)<360) 


%   g1=find(slice_ori(a)==0)
%   g2=find(slice_ori(a)==1)



% g1=find((out_ang_exL23(a,3)-out_ang_inL23(a,3))*69>0);
% g2=find((out_ang_exL23(a,3)-out_ang_inL23(a,3))*69<0);
% g1=find(od_out_iviv(:,2)<0.3)
% g2=find(od_out_iviv(:,2)>0.3)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2])
xticklabels({'45-135°','225-300'})
;xlabel('DIR');ylabel('L23 centroid EX-IN');
%% 
a=find(od_out_iviv(:,2)>0)  
%b=find(od_out_iviv(:,4)>135)
%[ia ib]=intersect(a,b)

fe=od_out_iviv(a,2)
centroid_plot(a,out_ang_exL23,out_ang_exL4,out_ang_exL5,out_ang_inL23,out_ang_inL4,out_ang_inL5,3,fe,{'DIR'})
%% 
corr_plot(od_out_iviv(:,2),od_out_iviv(:,5),(out_ang_exL23(:,3)-out_ang_exL23(:,1))*69,{'DIR','Horizontal centre of mass','h'});set(gca,'FontSize',12);
%% 
a=find(od_out_iviv(:,2)>0.2)  
b=find(od_out_iviv(:,5)>60 & od_out_iviv(:,5)<170) 
c=find(od_out_iviv(:,5)>225 & od_out_iviv(:,5)<300)
d=[b;c]
%c=find(pia_input<300)
[iaa ib]=intersect(a,d)
com=[];com=[L23fr(iaa,1) L4fr(iaa,1)  L23fr(iaa,2) L4fr(iaa,2) out_ang_exL23(iaa,5)  out_ang_exL4(iaa,5) out_ang_inL23(iaa,5)  out_ang_inL4(iaa,5) out_ang_exL23(iaa,3)  out_ang_exL4(iaa,3) out_ang_inL23(iaa,3)  out_ang_inL4(iaa,3) od_out_iviv(iaa,[1 2 3 4 5 6 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:20]);yticks([1:1:20]);

%% Display fraction for ex and in groups 
%group based on oripref
g1=find(od_out_iviv(:,8)<75);g2=find(od_out_iviv(:,8)>75);g3=[];
gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
gv(g1)=1;gv(g2)=2;gv(g3)=3;
gv(find(gv==0))=NaN;
%call function
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
g1=[];g2=[];g3=[];
%% Correlations Morphology and in vivo
com=[];com=[od_out_iviv morph_parameters(:,1:21) sum_densap' sum_densba' max_densba' max_densap' MLba_diff' MLap_diff']
correlation_matrix(com,0);title('PMorphology od out');
 %% %%Display desired correlations between input and iviv
corr_plot(com(:,4),com(:,19),com(:,12),{'L4ex input','Apical width/height','XSA'}); 
corr_plot(com(:,1),pia_input,pia_input,{'L4ex input','Apical width/height','XSA'}); 
%% 
corr_plot(com(find(com(:,7)<50),7),com(find(com(:,7)<50),13),pia_input(find(com(:,7)<50)),{'L4ex input','Apical width/height','XSA'}); 
%% 
find(com(:,12)>240)
corr_plot(com(find(com(:,12)>240),12),com(find(com(:,12)>240),4),pia_input(find(com(:,12)>240)),{'Maximal radial distance of apical tree from soma','Orientation preference','Pial depth'});  xlim([200 450]);ylim([0 200]);
%% 
discretize_plot(pia_input,3,od_out_iviv(:,4),1);xlabel('Pial depth bins');ylabel('Orientation preference');xticks([1:1:3]);title('in vivo in vitro cells');ylim([0 180])
%% 
corr_plot(com(:,23),com(:,1),pia_input,{'Number of branch points basal','Orientation selectivity index (gOSI)','Pial depth'}); xlim([0 22]);ylim([0 1])
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
com=[];com=[spon_out_iviv(:,:) out_ang_inL23]
correlation_matrix(com,0);title('ang od out');
 %% %%Display desired correlations between input and iviv
corr_plot(tot_inputL23(:,2),od_out_iviv(:,9),pia_input,{'L4ex input','Apical width/height','XSA'});
%% 
com=[];com=[sftf_out_iviv(:,[1:6]) tot_input]
correlation_matrix(com,0);title('ang od out');
 %% %%Display desired correlations between input and iviv
corr_plot(od_out_iviv(:,8),tot_input(:,2),od_out_iviv(:,8),{'L4ex input','Apical width/height','XSA'});
%% Plot Morphologies of cells all or defined subgroups directly from structure
id_m=nan_vector(morph_cells_id)
plot_morphologies(str,id_m,10,10);
%% Get all 147 morphtraces
for i=1:length(nan_vector)
    if ~isempty(str(nan_vector(i)).morph)==1 
    zz{:,i}=str(nan_vector(i)).morphtraces;
    else
        zz{:,i}=NaN;
    end
end

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
m_res=[];
for i=1:length(nan_vector)
    if ~isempty(str(nan_vector(i)).morph)==1 & str(nan_vector(i)).resp==1 & ~isempty(str(nan_vector(i)).OSI)==1
        m_res(i)=1;
     
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=find(m_res==1);
%% Plot cells sorted based on OSI
plot_morphologies_iviv(zz,morph_res_sub,6,6)
%% 
[isa isb]=sort(od_out_iviv(morph_res_sub,1));
%% 
plot_morphologies_iviv(zz,morph_res_sub(isb),6,6,m_flip_a)
%% 
plot_morphologies_iviv(zz,[123 135],6,6,m_flip_a)


%% Morphlogy to input correlations
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
%% Exication and Inhbition correlation of apical tree with map
figure;histogram([correlation_values_ap{morph_cells_id,1}],10);hold on;histogram([correlation_values_ap{morph_cells_id,2}],10);
%% 
for i=1:length(correlation_values_ap)
    if isempty(correlation_values_ap{i})==0
ap_correx(i)=[correlation_values_ap{i,1}];
ap_corrin(i)=[correlation_values_ap{i,2}];
    else
        ap_correx(i)=NaN;
        ap_corrin(i)=NaN;
    end
end
%% 
%% Display correlation between all input maps and in vivo
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) out_ang_exL23(:,5)  out_ang_exL4(:,5) out_ang_inL23(:,5)  out_ang_inL4(:,5) frh_lateral(:,1) frh_medial(:,2) od_out_iviv(:,[1 2 3 4 5 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);

%% Display correlation between all input maps and in vivo

com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) out_ang_exL23(:,5)  out_ang_exL4(:,5) out_ang_inL23(:,5)  out_ang_inL4(:,5) out_ang_exL23(:,3)  out_ang_exL4(:,3) out_ang_inL23(:,3)  out_ang_inL4(:,3) od_out_iviv(:,[1 2 3 7 8])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);

%% Display correlation between all input maps and in vivo

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 300]);imagesc(G(13:end,1:12));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:5]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','COM L2/3_{ex}','COM L4_{ex}','COM L2/3_{in}','COM L4_{in}'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','TW','Ca_{peak}'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% Subsample data OSI
a=find(od_out_iviv(:,1)>0.25)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) out_ang_exL23(a,5)  out_ang_exL4(a,5) out_ang_inL23(a,5)  out_ang_inL4(a,5) out_ang_exL23(a,3)  out_ang_exL4(a,3) out_ang_inL23(a,3)  out_ang_inL4(a,3) od_out_iviv(a,[1 2 3 4 5 6 7])]; 

correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
% xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});xtickangle(45)
% yticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});ytickangle(45)
%% ORI with cutoff
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 200]);imagesc(G(13:end,1:12));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:1]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','COM L2/3_{ex}','COM L4_{ex}','COM L2/3_{in}','COM L4_{in}'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'ORI'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);
%% Subsample data DSI
a=find(od_out_iviv(:,2)>0.3)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) out_ang_exL23(a,5)  out_ang_exL4(a,5) out_ang_inL23(a,5)  out_ang_inL4(a,5) out_ang_exL23(a,3)  out_ang_exL4(a,3) out_ang_inL23(a,3)  out_ang_inL4(a,3) od_out_iviv(a,[5])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
%% DIR with cutoff
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 200]);imagesc(G(13:end,1:12));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:1]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','COM L2/3_{ex}','COM L4_{ex}','COM L2/3_{in}','COM L4_{in}'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'DIR'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);


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
[idx_input, clustering_input, leafOrder] = hca([score_ex(:,pcs) score_in(:,pcs)],0,'ward',clu_num,pia_input,0,0.6);%call function for clustering
%[idx_input, clustering_input, leafOrder] = hca([data_w_input(:,pcs)],0,'ward',clu_num,pia_input,1);%call function for clustering
%% kmeans clustering
% idx_input=[];
% idx_input = kmeans([score_ex(:,pcs) score_in(:,pcs)],4)
%% Barplot difference of clusters
%Pial depth
[statsout] = barplot_sw(pia_input,idx_input,{'Clusters','Pial depth (µm)'});set(gca,'Ydir','reverse');xtickangle(45)
%% Orientation preference
[statsout] = barplot_sw(od_out_iviv(:,4),idx_input,{'Clusters','Orientation preference'});xtickangle(45);set(gca,'Ydir','reverse');set(gca,'FontSize',12)
%% Orientation selectivity
[statsout] = barplot_sw(od_out_iviv(:,8),idx_input,{'Clusters','Ca_{peak}'})
%% Angle centroid
[statsout] = barplot_sw(out_ang_inL23(:,3),idx_input,{'Clusters','Angle Centroid L23in'})
%% L4 ex fraction
[statsout] = barplot_sw(L4fr(:,1),idx_input,{'Clusters','L4 ex fraction'});set(gca,'FontSize',12)
%% L4 total input sum
[statsout] = barplot_sw(abs(tot_inputL4(:,1)),idx_input,{'Clusters','L4 total input'})
%% L23 in fraction
[statsout] = barplot_sw(L23fr(:,2),idx_input,{'Clusters','L23 in fraction'})
%% L23 in fraction
[statsout] = barplot_sw(tot_inputL23(:,2),idx_input,{'Clusters','L23 total input'})
%% 
[statsout] = barplot_sw(out_ang_inL23(:,5),idx_input,{'Clusters','Angle IN L23'})
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
dendroplot_SW(clustering_input,leafOrder,11,[score_ex(:,1:3) score_in(:,1:3)],{'PC1ex','PC2ex','PC3ex','PC1in','PC2in','PC3in'},pia_input)
%% Looking at input fractions
dendroplot_SW(clustering_input,leafOrder,11,[L23fr(:,1) L4fr(:,1)  L5fr(:,1) L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3ex','L4ex','L5ex','L2/3in','L4in','L5in'})
%% Looking at input fractions ex
dendroplot_SW(clustering_input,leafOrder,11,[L23fr(:,1) L4fr(:,1)  L5fr(:,1)],{'L2/3ex','L4ex','L5ex'})
%% Looking at input fractions in
dendroplot_SW(clustering_input,leafOrder,11,[L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3in','L4in','L5in'})
%% Plot clusters with values underneath it as heatmaps: SETUPS
dendroplot_SW(clustering_input,leafOrder,11,[setups],{'Setups'})
%% Plot clusters with values underneath it as heatmaps: SLICE ORIENTATION

dendroplot_SW(clustering_input,leafOrder,11,[slice_ori' score_in(:,3)],{'Slice Ori','PC2in'})
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
%% %% Morphology sholl analysis
close all;
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147);
%% 
ylabel('Number of dendritic crossing');xlabel('Distance from Soma (µm)');set(gca,'FontSize',12);
%% Correlations Morphology and in vivo
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[morph_parameters(:,2) nanmax(df,[],2) dis_s'  morph_parameters(:,4)  morph_parameters(:,5) max_s'  sum_densap' ...
    morph_parameters(:,12) nanmax(db,[],2) dis_s_ba'  morph_parameters(:,14)  morph_parameters(:,15) max_s_ba'  sum_densba'  od_out_iviv(:,[1 2 3 4 5 7 8])]
G=correlation_matrix(com,0);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 325]);imagesc(G(15:end,1:14));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:14]);yticks([1:1:7]);
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Max branch order','Peak number crossing','Total Density','Total Length','Max extent','Dis peak branch','Nr. branch points','Max branch order','Peak number crossing','Total Density'});xlabel('Feature');xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});ylabel('Feature');ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% 



%% 
%Ori and apical 
corr_plot(morph_parameters(:,3),od_out_iviv(:,4),pia_input,{'a','Apical width/height','Pial depth (µm)'});%xlim([-55 95]);ylim([0 180]);ylabel('Orientation preference','Color','k');xlabel('Angle centroid L2/3 IN','Color','b')
%corr_plot(L4fr(:,1),od_out_iviv(:,8),pia_input,{'L4ex input','Apical width/height','XSA'});
%% 
plot_morphologies_iviv(zz,[123 135],6,6,m_flip_a);
%% 
corr_plot(morph_parameters(id_ori,3),out_ang_inL23(id_ori,5),pia_input(id_ori),{'a','Apical width/height','Pial depth (µm)'})
%% OSI and morph apical
corr_plot(morph_parameters(morph_res_sub,4),od_out_iviv(morph_res_sub,1),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});xlim([0 25]);ylim([0 1]);ylabel('gOSI','Color','k');xlabel('Number of branch points apical','Color','k');set(gca,'FontSize',12)
yticks(0:0.25:1)
%% OSI and morph basal
corr_plot(morph_parameters(morph_res_sub,14),od_out_iviv(morph_res_sub,1),pia_input(morph_res_sub),{'a','Apical width/height','Pial depth (µm)'});xlim([0 30]);ylim([0 1]);ylabel('gOSI','Color','k');xlabel('Number of branch points basal','Color','k');set(gca,'FontSize',12)
yticks(0:0.25:1)
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
a=99;
ml=15;
  tmp=zz{1,a}
  figure;set(gcf,'color','w')
  m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
         m.EdgeColor = [0 0 0]
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
 
 