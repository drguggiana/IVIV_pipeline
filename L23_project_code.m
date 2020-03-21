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
    %get the layers
   % lyr = invitro_struct(cells).layers;
    %get the overlap map
    ov_map(:,:,cells) = squeeze(sum(cat(3,bin_excm,bin_inh),3));  
    ov_map_n(:,:,cells)=ov_map(:,:,cells)./max(max(ov_map(:,:,cells)));
    %Get the ex and in total per map
    ex_tot(:,cells)=str(nan_vector(cells)).excinhTotal(1);
    in_tot(:,cells)=str(nan_vector(cells)).excinhTotal(2);
    %layers(:,:,cells)=str(nan_vector(i)).layers;
end
%% 
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
set(gca,'Ydir','reverse')
%histogram(pia_input,'FaceColor','k','FaceAlpha',0.1)
set(gca,'FontSize',14);
%set(gca,'FontWeight','bold')
%% Setup A and Setup B
setups=[zeros(47,1);ones(100,1)]

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
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv pialD layer_assign] = iviv_profiles(nan_vector(incl_idx:end),str);
frac_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv];
abs_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv];
layer_assign=[zeros(length(nan_vector(incl_idx:end)),2) layer_assign'];
%% Calculate the fraction from the overlap for verti and hori
[frac_ovh abs_ovh frac_ovv abs_ovv] = iviv_profiles_ov(ov_map);
%% Calculate fraction using the overlap binary maps
[frac_bh abs_bh frac_bv abs_bv] = iviv_profiles_ov(ov_map_bin);
%% Calculate difference between ex and in (ex-in) for L23, L4,  L5 ONE APPROACH
% for i=1:length(nan_vector(incl_idx:end))
% diffL23(i)=sum(sum(diff_map(3:5,:,i),2))/length(nonzeros(diff_map(3:5,:,i)));
% diffL4(i)=sum(sum(diff_map(6:7,:,i),2))/length(nonzeros(diff_map(6:7,:,i)));
% diffL5(i)=sum(sum(diff_map(8:11,:,i),2))/length(nonzeros(diff_map(8:11,:,i)));
% end
%% Calculate difference between ex and in (ex-in) for L23, L4,  L5 using fractions USING THIS FOR NOW
frac_diffv=frac_exv_m-frac_inv;
frac_diffh=frac_exh-frac_inh;
diffL23fr=nanmean(frac_diffv(:,3:5),2);
diffL4fr=nanmean(frac_diffv(:,6:7),2);
diffL5fr=nanmean(frac_diffv(:,8:10),2);
%% Calculate VERTICAL ex and in fraction for L23, L4,  L5
L23fr=[nanmean(frac_exv_m(:,3:5),2) nanmean(frac_inv(:,3:5),2)];
L4fr=[nanmean(frac_exv_m(:,6:7),2) nanmean(frac_inv(:,6:7),2)];
L5fr=[nanmean(frac_exv_m(:,8:10),2) nanmean(frac_inv(:,9:11),2)];
L23frov=[nanmean(frac_ovv(:,3:5),2)];
L4frov=[nanmean(frac_ovv(:,6:7),2)];
L5frov=[nanmean(frac_ovv(:,8:10),2)];
%% Calculate fraction using the layer assignement
for i=1:length(frac_exv_m)
    L4fr_l(i,:)=[nanmean(frac_exv_m(i,find(layer_assign(i,:)==3))) nanmean(frac_inv(i,find(layer_assign(i,:)==3)))]
    
end
%% 
figure;set(gcf,'color','w');imagesc(layer_assign')
hold on;line([47 47], [1 16],'Color','k','LineStyle','--');set(gca,'FontSize',10);
ylim([1 16]);yticks([1:1:16]);ylabel('Rows');xlabel('Cells');%c=colorbar;
%% Calculate difference between ex and in (ex-in) medial and lateral
frh_medial=[nanmean(frac_exh(:,1:8),2) nanmean(frac_inh(:,1:8),2)];
frh_lateral=[nanmean(frac_exh(:,9:end),2) nanmean(frac_inh(:,9:end),2)];
frh_diff_medial=nanmean(frac_diffh(:,1:8),2);
frh_diff_lateral=nanmean(frac_diffh(:,9:end),2);
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
for i=1:length(nan_vector(incl_idx:end))
tmp1=ex_map(3:5,:,i);
for k=1:16
tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
end
tmp2=in_map(3:5,:,i);
for k=1:16
tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
end
tmp=find(tL23e>0);
ex_spanhL23(i)=tmp(end)-tmp(1);
tmp=[];
tmp=find(tL23i>0);
in_spanhL23(i)=tmp(end)-tmp(1);
tmp=[];
tL23i=[];
tL23e=[];
tmp1=[];tmp2=[];
end
for i=1:length(nan_vector(incl_idx:end))
tmp1=ex_map(6:7,:,i);
for k=1:16
tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
end
tmp2=in_map(6:7,:,i);
for k=1:16
tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
end
tmp=find(tL23e>0);
if isempty(tmp)==0
ex_spanhL4(i)=tmp(end)-tmp(1);
else
ex_spanhL4(i)=0;
end
tmp=[];
tmp=find(tL23i>0);
if isempty(tmp)==0
in_spanhL4(i)=tmp(end)-tmp(1);
else
in_spanhL4(i)=0;
end
tmp=[];
tL23i=[];
tL23e=[];
end
%L5
for i=1:length(nan_vector(incl_idx:end))
tmp1=ex_map(8:10,:,i);
for k=1:16
tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
end
tmp2=in_map(8:10,:,i);
for k=1:16
tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
end
tmp=find(tL23e>0);
if isempty(tmp)==0
ex_spanhL5(i)=tmp(end)-tmp(1);
else
ex_spanhL5(i)=0;
end
tmp=[];
tmp=find(tL23i>0);
if isempty(tmp)==0
in_spanhL5(i)=tmp(end)-tmp(1);
else
in_spanhL5(i)=0;
end
tmp=[];
tL23i=[];
tL23e=[];
end
%% total sum absolute for ex and inh in pA
for i=1:length(nan_vector(incl_idx:end))
%whole map
tmp1=ex_map_raw(:,:,i);
tmp2=in_map_raw(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map_raw(3:5,:,i);
tmp2=in_map_raw(3:5,:,i);
tot_inputL23(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map_raw(6:7,:,i);
tmp2=in_map_raw(6:7,:,i);
tot_inputL4(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
end
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,[]);
%% Display fraction for ex and in as well as overlap for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_ovv*-1,frac_ovh*-1,[]);
%% Display fraction for ex and in as well as overlap bin for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_bv*-1,frac_bh*-1,[]);
%% Display both maximum horizontal span and diff betwen ex and in ONE APPROACH
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 200]);
%Diff between ex and in shwon with histogram and CDF for L23 L4 and L5
left_color =[0 0 0];right_color =[0 0 0]
subplot(1,2,1);
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);hold on;
h4 = histogram(diffL23fr,8);h4.BinWidth = 0.025;h4.EdgeColor = 'k';h4.FaceColor = 'm';hold on;h4 = histogram(diffL4fr,8);h4.BinWidth =  0.025;h4.EdgeColor = 'k';h4.FaceColor = 'g'
hold on;h4 = histogram(nonzeros(diffL5fr),8);h4.BinWidth =  0.025;h4.FaceColor = [0.5 0.5 0.5];ylabel('Cell counts');box off;
%yyaxis right;p1=cdfplot(diffL23fr);hold on;p2=cdfplot(diffL4fr);hold on;p3=cdfplot(nonzeros(diffL5fr));grid off; title('');
%ylabel('Cumulative');xlabel('Ex-In');p1.Color='m';p2.Color=[0 1 0];p3.Color=[0.5 0.5 0.5];p1.LineStyle='--';p2.LineStyle='--';p3.LineStyle='--';
xlim([-0.3 0.3]);xlabel('Ex-In');hold on;set(gca,'FontSize',10);hold on; title('Vertical');
subplot(1,2,2);
h4 = histogram(ex_spanhL23-in_spanhL23,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'm';
hold on;h4 = histogram(ex_spanhL4-in_spanhL4,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'g';
hold on;h4 = histogram(ex_spanhL5-in_spanhL5,8);h4.BinWidth = 1;box off;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];xlabel('\Delta Span');ylabel('Cell counts');
%yyaxis right;p1=cdfplot(ex_spanhL23-in_spanhL23);hold on;p2=cdfplot(ex_spanhL4-in_spanhL4);hold on;p3=cdfplot(ex_spanhL5-in_spanhL5);grid off; title('');
%ylabel('Cumulative');xlabel('Ex-In');p1.Color='m';p2.Color=[0 1 0];p3.Color=[0.5 0.5 0.5];legend('L23', 'L4','L5');legend boxoff; set(gca,'FontSize',10);
xlim([-10 10]);hold on; title('Horizontal');%p1.LineStyle='--';p2.LineStyle='--';p3.LineStyle='--';
legend('L23', 'L4','L5');legend boxoff; set(gca,'FontSize',10);
%% DISPLAY SPAN, SECOND APPROACH
display_sortfr(L23fr,1,'L2/3')
display_sortfr(L4fr,1,'L4')
display_sortfr(L5fr,1,'L5')
%% 
spanhL23=[ex_spanhL23;in_spanhL23]'
spanhL4=[ex_spanhL4;in_spanhL4]'
spanhL5=[ex_spanhL5;in_spanhL5]'
 
display_sortfr(spanhL23,2,'L2/3');
display_sortfr(spanhL4,2,'L4');
display_sortfr(spanhL5,2,'L5');

%% Alternative displaying EX and IN 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 800, 200]);
subplot(1,2,1);
binRange = -0.3:0.075:0.3;
hcx = histcounts(diffL23fr*-1,[binRange Inf]);
hcy = histcounts(diffL4fr*-1,[binRange Inf]);
hcz = histcounts(diffL5fr*-1,[binRange Inf]);
b1=bar(binRange,[hcx;hcy;hcz]');box off,
b1(1).FaceColor=[0 1 0];
b1(2).FaceColor='m'
b1(3).FaceColor=[0.8 0.8 0.8];
xlim([-0.35 0.35]);xlabel('EX - IN');hold on;set(gca,'FontSize',10);hold on; title('Vertical');ylabel('Cell counts');
ylim([0 100])
text(-0.15,100,'EX','Color','r');
text(0.15,100,'IN','Color','b');

subplot(1,2,2);
binRange = -8:2:10;
hcx = histcounts((ex_spanhL23-in_spanhL23)*-1,[binRange Inf]);
hcy = histcounts((ex_spanhL4-in_spanhL4)*-1,[binRange Inf]);
hcz = histcounts((ex_spanhL5-in_spanhL5)*-1,[binRange Inf]);
b1=bar(binRange,[hcx;hcy;hcz]');box off
b1(1).FaceColor=[0 1 0];
b1(2).FaceColor='m'
b1(3).FaceColor=[0.8 0.8 0.8];
xlim([-10.5 10.5]);xlabel('EX - IN');hold on;set(gca,'FontSize',10);hold on; title('Horizontal');ylabel('Cell counts');
ylim([0 80]);
text(-6,80,'EX','Color','r');
text(6,80,'IN','Color','b');
legend('L23', 'L4','L5');legend boxoff; set(gca,'FontSize',10);

%% Display both maximum horizontal span and diff betwen ex and in COLUMN APPROACH
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 450, 500]);
subplot(3,2,1);
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);hold on;
h4 = histogram(diffL23fr,8);h4.BinWidth = 0.025;h4.EdgeColor = 'k';h4.FaceColor = 'm';hold on;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');set(gca,'FontSize',10);
title('Vertical');ylabel('Cell counts');hold on;set(gca,'FontSize',10);
xlim([-0.3 0.3]);text(-0.2,60,'IN','Color','b');text(0.2,60,'EX','Color','r');text(-0.2,50,'L2/3','Color','k');
subplot(3,2,3);
h4 = histogram(diffL4fr,8);h4.BinWidth =  0.025;h4.EdgeColor = 'k';h4.FaceColor = 'g'
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');set(gca,'FontSize',10);ylabel('Cell counts');hold on;set(gca,'FontSize',10);box off;
text(-0.2,50,'L4','Color','k');
xlim([-0.3 0.3]);
subplot(3,2,5);
h4 = histogram(nonzeros(diffL5fr),8);h4.BinWidth =  0.025;h4.FaceColor = [0.5 0.5 0.5];ylabel('Cell counts');box off;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');
xlim([-0.3 0.3]);xlabel('Ex-In');hold on;set(gca,'FontSize',10);
text(-0.2,50,'L5','Color','k');

subplot(3,2,2);
h4 = histogram(ex_spanhL23-in_spanhL23,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'm';
xlim([-10 10]);hold on; title('Horizontal');box off;
hold on;line([0 0], [1 60],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7,60,'IN','Color','b');text(7,60,'EX','Color','r');
text(-7,50,'L2/3','Color','k');
subplot(3,2,4);
hold on;h4 = histogram(ex_spanhL4-in_spanhL4,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'g';
xlim([-10 10]);hold on;line([0 0], [1 60],'Color','k','LineStyle','--');set(gca,'FontSize',10);
text(-7,50,'L4','Color','k');
subplot(3,2,6);
hold on;h4 = histogram(ex_spanhL5-in_spanhL5,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];;
xlim([-10 10]);hold on;line([0 0], [1 60],'Color','k','LineStyle','--');set(gca,'FontSize',10);
xlabel('Ex-In');hold on;set(gca,'FontSize',10);
text(-7,50,'L5','Color','k');
%% 
subplot(3,2,2);
h4 = histogram(ex_spanhL23-in_spanhL23,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'm';
hold on;h4 = histogram(ex_spanhL4-in_spanhL4,8);h4.BinWidth = 1;h4.EdgeColor = 'k';h4.FaceColor = 'g';
hold on;h4 = histogram(ex_spanhL5-in_spanhL5,8);h4.BinWidth = 1;box off;h4.EdgeColor = 'k';h4.FaceColor = [0.5 0.5 0.5];xlabel('\Delta Span');ylabel('Cell counts');
%yyaxis right;p1=cdfplot(ex_spanhL23-in_spanhL23);hold on;p2=cdfplot(ex_spanhL4-in_spanhL4);hold on;p3=cdfplot(ex_spanhL5-in_spanhL5);grid off; title('');
%ylabel('Cumulative');xlabel('Ex-In');p1.Color='m';p2.Color=[0 1 0];p3.Color=[0.5 0.5 0.5];legend('L23', 'L4','L5');legend boxoff; set(gca,'FontSize',10);
xlim([-10 10]);hold on; title('Horizontal');%p1.LineStyle='--';p2.LineStyle='--';p3.LineStyle='--';
legend('L23', 'L4','L5');legend boxoff; set(gca,'FontSize',10);
%% binary overview
fig1=figure;set(gcf,'color','w');hold on;
for i=1:size(frac_bv,1)
exp=plot(frac_bv(i,1:16)',1:16,'-m');set(gca,'Ydir','reverse');ylabel('Vertical input');xlabel('Fraction of overlap EX/IN');box off;
exp.Color(4) = 0.2;
end
hold on;mexp=errorbar(nanmean(frac_bv(:,1:16)),1:16,nanstd(frac_bv(:,1:16))/sqrt(size(frac_bv,1)),'horizontal','k');set(gca,'Ydir','reverse');
hold on;line([0 0.5], [3 3],'Color','k','LineStyle','--');hold on;line([0 0.5], [6 6],'Color','k','LineStyle','--');
hold on;line([0 0.5], [8 8],'Color','k','LineStyle','--');hold on;line([0 0.5], [11 11],'Color','k','LineStyle','--');hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});
set(gca,'FontSize',12);xlim([0 0.5])
fig1=figure;set(gcf,'color','w');hold on;
for i=1:size(frac_bv,1)
exp=plot(frac_bh(i,1:16)','-m');xlabel('Horizontal input');ylabel('Fraction of overlap EX/IN');box off;
exp.Color(4) = 0.2;
end
hold on;mexp=errorbar(nanmean(frac_bh(:,1:16)),nanstd(frac_bh(:,1:16))/sqrt(size(frac_bh,1)),'horizontal','k');
hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8 8], [0 0.5],'Color','k','LineStyle','--');set(gca,'FontSize',10);set(gca,'FontSize',12);ylim([0 0.5])
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
%group based on pial depth
g1=find(pia_input>=220);g2=find(pia_input<220);g3=[];
gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
gv(g1)=1;gv(g2)=2;gv(g3)=3;
%call function
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
g1=[];g2=[];g3=[];
%% %%Display desired correlations between PCs and actual data 
%EX
corr_plot(L23fr(:,1),score_ex(:,1),pia_input,{'PC1ex','L23ex','Pial depth'});
corr_plot(L4fr(:,1),score_ex(:,1),pia_input,{'PC1ex','L4ex','Pial depth'});
corr_plot(nonzeros(L5fr(:,1)),score_ex(find(nonzeros(L5fr(:,1))),1),pia_input(find(nonzeros(L5fr(:,1)))),{'PC1ex','L5ex','Pial depth'});
%IN
corr_plot(L23fr(:,2),score_in(:,1),pia_input,{'PC1in','L23in','Pial depth'});
corr_plot(L4fr(:,2),score_in(:,1),pia_input,{'PC1in','L4in','Pial depth'});
corr_plot(nonzeros(L5fr(:,2)),score_in(find(nonzeros(L5fr(:,2))),1),pia_input(find(nonzeros(L5fr(:,2)))),{'PC1in','L5in','Pial depth'});
%% %%Display desired correlations between pial depth and L23/L4
%EX
corr_plot(L23fr(:,1),pia_input,[],{'L23ex input','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(L4fr(:,1),pia_input,[],{'L4ex input','Pial depth',;});set(gca,'Ydir','reverse');
corr_plot(tot_inputL4(:,2),pia_input,[],{'L4 total input ex','Pial depth',;});set(gca,'Ydir','reverse');
corr_plot(tot_inputL23(:,2),pia_input,[],{'L23 total input in','Pial depth',;});set(gca,'Ydir','reverse');
%IN
corr_plot(L23fr(:,2),pia_input,[],{'L23in input','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(L4fr(:,2),pia_input,[],{'L4in input','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(nonzeros(L5fr(:,1)),pia_input(find(nonzeros(L5fr(:,1)))),pia_input(find(nonzeros(L5fr(:,1)))),{'PC1in','L5in','Pial depth'});
%% Pial depth correlation plot with fraction 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 600, 200]);
%Diff between ex and in shwon with histogram and CDF for L23 L4 and L5
subplot(1,3,1);
par1=L4fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))

[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;xlabel('L4 fraction');ylabel('Pial depth (µm)'); 
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'r-');set(gca,'box','off');set(gcf,'color','w');axis square;   set(gca,'Ydir','reverse');
hold on;
par1=L4fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;xlabel('L4 fraction');ylabel('Pial depth (µm)'); set(gca,'FontSize',10)
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'b-');set(gca,'box','off');set(gcf,'color','w');axis square;   set(gca,'Ydir','reverse');
    %legend('Ex','In')
  
    
subplot(1,3,2);
par1=L23fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))

[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;xlabel('L23 fraction');ylabel('Pial depth (µm)'); set(gca,'FontSize',10)
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
  P = polyfit(par1,par2,1);
    yfit = P(1)*par1+P(2);
    hold on;
    plot(par1,yfit,'r-');set(gca,'box','off');set(gcf,'color','w');axis square;   set(gca,'Ydir','reverse');
hold on;
par1=L23fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;xlabel('L23 fraction');ylabel('Pial depth (µm)'); 
 if P(2)<0.05 & P(2)>0.01
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.05'])
     
 elseif P(2)<0.01 & P(2)>0.001
     title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.01'])
 elseif P(2)<0.001
    title(['r= ' mat2str(round(R(2),2)) ' ' 'p<0.001'])
 else
     title(['r= ' mat2str(round(R(2),2)) ' ' 'n.s'])
 end
 
subplot(1,3,3);
par1=L5fr(:,1);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))

[R P]=corrcoef(par1,par2);s=scatter(par1,par2,5,'o','MarkerEdgeColor','r','MarkerFaceColor','r');box off;xlabel('L5 fraction');ylabel('Pial depth (µm)'); set(gca,'FontSize',10)
set(gca,'Ydir','reverse');
hold on;
par1=L5fr(:,2);
par2=pia_input((find(par1>0)));
par1=par1((find(par1>0)))
[R P]=corrcoef(par1,par2);scatter(par1,par2,5,'o','MarkerEdgeColor','b','MarkerFaceColor','b');box off;xlabel('L5 fraction');ylabel('Pial depth (µm)'); 
 set(gca,'Ydir','reverse');xlim([0 0.4]);axis square
 
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
%EX ang centroid
[out_ang_ex] = centroid_map(ex_map(:,:,:),somax,pia_input,[1:cells],0);
[out_ang_exL23] = centroid_map(ex_map(3:5,:,:),somax,pia_input,[1:cells],2);
[out_ang_exL4] = centroid_map(ex_map(6:7,:,:),somax,pia_input,[1:cells],5);
%difference map
[out_ang_diff] = centroid_map(diff_map(3:5,:,:),somax,pia_input,[1:cells],0);
%fake map
[out_ang_fake] = centroid_map(map_fake(3:5,:,:),somax,pia_input,[1:cells],2);
%% Aligned maps
for i=1:length(nan_vector)
align_in(:,:,i)=reshape(aligned_maps_in(i,:),22,16);
end
%% 

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
%% DIstributions of actual centroid
ang2=out_ang_inL23;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 900, 200]);
subplot(1,4,1);
h4 = histogram(ang2(:,4)*69-ang2(:,2)*69,4);h4.BinWidth = 0.5*69;box off;h4.EdgeColor = 'k';h4.FaceColor = 'w';xlabel('\Delta Centroid Y');ylabel('counts')
subplot(1,4,2);
scatter(ang2(:,4)*69,ang2(:,2)*69,'ko');set(gcf,'color','w');[R P]=corrcoef(ang2(:,4),ang2(:,2),'row','complete');xlabel('Vertical Centroid (\mum)');ylabel('Soma location y (\mum)')
subplot(1,4,3);
h4 = histogram(ang2(:,3)*69-ang2(:,1)*69,4);h4.BinWidth = 0.5*69;box off;h4.EdgeColor = 'k';h4.FaceColor = 'w';xlabel('\Delta Centroid X');ylabel('counts')
subplot(1,4,4);
scatter(ang2(:,3)*69,ang2(:,2)*69,'ko');set(gcf,'color','w');xlabel('Horizontal Centroid (\mum)');ylabel('Soma location x (\mum)');
%% Relation betwen ex and in centroid
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
%% Display correlation between all PCs and pial depth Multiple comparison 
com=[score_ex(:,1:3) score_in(:,1:3) L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov L5frov out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_diff(:,5)]; 
correlation_matrix(com,1);title('Vertical');
com=[score_ex(:,1:3) score_in(:,1:3) frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)']; 
correlation_matrix(com,1);title('Horizontal');
%% Display correlations scores and angle towards centroid
com=[score_ex(:,1:3) score_in(:,1:3) out_ang_inL23 pia_input]; 
correlation_matrix(com,1);title('PC scores and angle to centroid');
%% %%Display desired correlations PC2in and angle and pia input
corr_plot(out_ang_inL23(:,3),score_in(:,2),diffL4',{'Angle centroid L23 in','PCin2','L4ex-in'});set(gca,'Ydir','reverse');
%% Desired correlation between angle and inhbition 
corr_plot(out_ang_inL23(:,5),score_in(:,2),pia_input,{'Angle centroid L2/3 IN','PC2_{in}','Pial depth (µm)'});ylabel('PC2_{in}','Color','b');xlabel('Angle centroid L2/3 IN','Color','b')
%% Desired correlation between angle and pial depth 
corr_plot(out_ang_inL23(:,4),score_in(:,2),pia_input,{'Y centroid in L23in','PCin2','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(out_ang_inL23(:,5),pia_input,pia_input,{'Angle centroid L23 in','Pial depth','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(out_ang_exL4(:,5),pia_input,pia_input,{'Angle centroid L23 ex','Pial depth','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(L23fr(:,2),pia_input,out_ang_inL23(:,5),{'L23in fraction','Pial depth','angle L23in'});set(gca,'Ydir','reverse');
%% FAKE MAPS CENTROID
corr_plot(out_ang_fake(:,5),score_in(:,2),pia_input,{'Angle centroid L23 fake','PCin2','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(out_ang_fake(id_ori,5),score_in(id_ori,2),pia_input(id_ori),{'Angle centroid L23 fake','PCin2','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(out_ang_fake(:,5),pia_input,pia_input,{'Angle centroid L23 fake','Pial depth','Pial depth'});set(gca,'Ydir','reverse');
%% 

corr_plot(out_ang_inL23(id_ori,5),score_in(id_ori,2),pia_input(id_ori),{'Angle centroid L2/3 IN','PC2_{in}','Pial depth (µm)'});set(gca,'Ydir','reverse');

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
correlation_matrix(com,0);title('Morphology vs PC input scores');
%% %% Morphology parameters vs real parameters vertical
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov tot_inputL23 tot_inputL4 morph_parameters(:,1:21) sum_densap' sum_densba']; 
correlation_matrix(com,1);title('Morphology vs Vertical input');
%% %% Morphology parameters vs real parameters angle
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_diff(:,5) morph_parameters(:,1:21) sum_densap' sum_densba']; 
correlation_matrix(com,1);title('Morphology vs Vertical input');
%% %% Morphology parameters vs real parameters horizontal
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff' sum_densap' sum_densba']; 
correlation_matrix(com,0);title('Morphology vs Horizontal input');
%% %%Display desired correlations between pial depth and L23/L4
corr_plot(L4fr(:,1),morph_parameters(:,9),pia_input,{'L4ex input','Apical width/height','Pial depth'});
%% %%Display desired correlations between pial depth and L23/L4
corr_plot(diffL23,morph_parameters(:,13),pia_input,{'L4 ex/in ov','Max Apical','Pial depth'});
%% Morphology and input map correlations
corr_plot(out_ang_inL23(:,5),morph_parameters(:,1),pia_input,{'Angle centroid L23','Apical vertical span','Pial depth'});set(gca,'Ydir','reverse');
%% Using the Morph PCs for correlation
scores_m=ones(148,3)*NaN;
scores_m(morph_cells_id,:)=score_morph(:,1:3);
%% %% Morphology PCs vs PCs input
com=[];com=[score_ex(:,1:3) score_in(:,1:3) scores_m]; 
correlation_matrix(com,1);title('PCscors input vs Morph scores');
%% %% Morphology PCs vs PCs input
com=[];com=[scores_m L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral]; 
correlation_matrix(com,1);title('PCscors input vs Morph scores');
%% 
% %% Display fraction for ex and in as well as diff for all 16 rows and columns based on TMD groups
% load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\com_189_withTMDclusteridx.mat');
% [val idxt idxm] = intersect(cellID_str,com_TMDidx(:,28))
% tmd_gr=com_TMDidx(idxm,27);
% g1=find(tmd_gr==1);g2=find(tmd_gr==2);g3=[];
% gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
% gv(g1)=1;gv(g2)=2;gv(g3)=3;
% %call function
% [stats_g] = display_inputs([frac_exv_m(idxt,:) frac_inv(idxt,:)],[frac_exh(idxt,:) frac_inh(idxt,:)],frac_diffv(idxt,:),frac_diffh(idxt,:),gv);
% g1=[];g2=[];g3=[];

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


%%  IN VIVO IN VITRO
%% Read out paramteres from iviv structure
[od_out_iviv spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(str,nan_vector);

%% Get IDS for ORI
id_ori=find(~isnan(od_out_iviv(:,4)));
%% Correlations PCs input and OD out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) od_out_iviv]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations PCs input and sftft out
com=[];com=[score_ex(:,1:3) score_in(:,1:3) sftf_out_iviv(:,:)]; 
correlation_matrix(com,0);title('PC inputs ODout iviv');
%% Correlations Input vertical vs ODout
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23fr diffL4fr L23frov L4frov L4fr(:,1)-L23fr(:,2) tot_input tot_inputL23 tot_inputL4 od_out_iviv(:,1:8)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Correlations Input vertical vs ODout
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) tot_input(:,1) od_out_iviv]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');
%% Correlations Input vertical vs ODout
com=[];com=[out_ang_ex(:,5) out_ang_in(:,5) out_ang_inL23(:,5) out_ang_exL23(:,5) out_ang_inL4(:,5) out_ang_exL4(:,5) out_ang_exL4(:,5)-out_ang_inL23(:,5) out_ang_diff(:,5) sftf_out_pref_iviv(:,1:3)]; 
correlation_matrix(com,0);title('Input vertical ODout iviv');

 %% %%Display desired correlations between pial depth and L23/L4
corr_plot(tot_input(:,1),sftf_out_iviv(:,8),od_out_iviv(:,end),{'L4ex input','Apical width/height','XSA'});
%%  %% %%Display desired correlations between pial depth and L23/L4
corr_plot(score_in(:,2),od_out_iviv(:,4),pia_input,{'PC2in','Orientation preference','Pial depth'});
%%  %% %%Display desired correlations between pial depth and L23/L4
corr_plot(out_ang_inL23(:,10),od_out_iviv(:,4),pia_input,{'Angle centroid L23in','Orientation preference','Pial depth'});%xlim([-60 120]) 
%%  %% %%Display desired correlations between pial depth and L23/L4
corr_plot(out_ang_exL4(:,5),od_out_iviv(:,4),L4fr(:,1),{'Angle centroid L4ex','Orientation preference','Pial depth'});;xlim([45 90]) 
%% Display correlation with fake map
corr_plot(out_ang_fake(:,5),od_out_iviv(:,4),pia_input,{'Angle centroid L4ex','Orientation preference','Pial depth'});
%% 
corr_plot(out_ang_exL23(:,4),od_out_iviv(:,4),pia_input,{'Angle centroid L4ex','Orientation preference','Pial depth'});
%% 
corr_plot(out_ang_inL23(:,2)-out_ang_exL4(:,4),od_out_iviv(:,4),pia_input,{'Angle centroid L4ex','Orientation preference','Pial depth'});
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
correlation_matrix(com,1);title('Input horizontal spon_out iviv');
%% Correlations PCs input and sftf_out_iviv
com=[];com=[score_ex(:,1:3) score_in(:,1:3) sftf_out_iviv(:,:)]; 
correlation_matrix(com,0);title('PC inputs sftf_out iviv');
%% Correlations PCs input and sftf_out_iviv
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23' diffL4' L23frov L4frov sftf_out_iviv(:,[1:6])]; 
correlation_matrix(com,0);title('PC inputs sftf_out iviv');
%% Correlations PCs input and sftf_out_iviv
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' sftf_out_iviv]; 
correlation_matrix(com,1);title('PC inputs sftf_out iviv');
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
%% 

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 500, 325]);imagesc(G(11:end,1:10));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:10]);yticks([1:1:7]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}'});xlabel('Feature');xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});ylabel('Feature');ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)

% xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});xtickangle(45)
% yticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','lat_{ex}','med_{in}','OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});ytickangle(45)
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

%%  morphology
[statsout] = barplot_sw(MLba_diff',idx_input,{'Clusters','L23 total input'})
%% 
[statsout] = barplot_sw(ap_corrin',idx_input,{'Clusters','Correlation morpho input'})
%% 
[statsout] = barplot_sw(morph_parameters(id_ori,3),idx_input(id_ori),{'Clusters','Correlation morpho input'})
%% Plot clusters with values underneath it as heatmaps
dendroplot(clustering_input,leafOrder,11,[score_ex(:,1:3) score_in(:,1:3)],{'PC1ex','PC2ex','PC3ex','PC1in','PC2in','PC3in'},pia_input)
%% Looking at input fractions
dendroplot(clustering_input,leafOrder,11,[L23fr(:,1) L4fr(:,1)  L5fr(:,1) L23fr(:,2) L4fr(:,2)  L5fr(:,2)],{'L2/3ex','L4ex','L5ex','L2/3in','L4in','L5in'})
%% Plot clusters with values underneath it as heatmaps
dendroplot(clustering_input,leafOrder,11,[setups],{'Setups'})


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
a=94;
ml=15;
  tmp=zz{1,a}
  figure;set(gcf,'color','w')
  m=plot_tree(tmp{1,1},[1 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
         m.EdgeColor = [0 0 0]
       m1=plot_tree(tmp{1,2},[0 0 0],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = 'm'
       m2=plot_tree(tmp{1,3},[0 0 1],[0 tmp{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = 'g'
set(gca,'Ydir','reverse');
 

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
 
 