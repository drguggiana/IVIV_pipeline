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
%% Display variance explained for ex and in maps
var_exp(explained_in,explained_ex,{'Inhibition','Excitation'});  
%% %% Display variance explained combined
var_exp(explained_com,[],[]); 
%% %% Display variance explained diff maps
var_exp(explained_diff,[],[]); 
%% %% Display variance explained basal and apical
var_exp(explained_basal,explained_apical,{'Basal','Apical'}); 
%% Display coefficent of PCs ALIGNED in and ex
coeff_display(coeff_ex,coeff_in,bin_num,hbin_num);
%% Display coefficent of PCs ALIGNED combined
coeff_display(coeff_com(1:352,:),coeff_com(353:end,:),bin_num,hbin_num);
%% Display coefficent of PCs ALIGNED diff
coeff_display(coeff_diff,[],bin_num,hbin_num);
%% Display correlation between all PCs and pial depth Multiple comparison 
com=[score_ex(:,1:3) score_in(:,1:3) score_com(:,1:3) score_diff(:,1:3) pia_input]; 
correlation_matrix(com,1);
%% %%Display desired correlations between PCs/pia
corr_plot(score_ex(:,1),score_com(:,1),pia_input,{'PC1ex','PC1com','Pial depth'});
%% Load fraction and absolute input of maps for ex and inh
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv pialD layer_assign] = iviv_profiles(nan_vector(incl_idx:end),str);
frac_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv];
abs_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv];
%% Calculate the fraction from the overlap for verti and hori
[frac_ovh abs_ovh frac_ovv abs_ovv] = iviv_profiles_ov(ov_map);
%% Calculate difference between ex and in (ex-in) for L23, L4,  L5
for i=1:length(nan_vector(incl_idx:end))
diffL23(i)=sum(sum(diff_map(3:5,:,i),2))/length(nonzeros(diff_map(3:5,:,i)));
diffL4(i)=sum(sum(diff_map(6:7,:,i),2))/length(nonzeros(diff_map(6:7,:,i)));
diffL5(i)=sum(sum(diff_map(8:11,:,i),2))/length(nonzeros(diff_map(8:11,:,i)));
end
%% Calculate difference between ex and in (ex-in) for L23, L4,  L5 using fractions
frac_diffv=frac_exv_m-frac_inv;
frac_diffh=frac_exh-frac_inh;
diffL23fr=nanmean(frac_diffv(:,3:5),2);
diffL4fr=nanmean(frac_diffv(:,6:7),2);
diffL5fr=nanmean(frac_diffv(:,8:10),2);
%% Calculate VERTICAL ex and in fraction for L23, L4,  L5
L23fr=[nanmean(frac_exv_m(:,3:5),2) nanmean(frac_inv(:,3:5),2)];
L4fr=[nanmean(frac_exv_m(:,6:7),2) nanmean(frac_inv(:,6:7),2)];
L5fr=[nanmean(frac_exv_m(:,8:10),2) nanmean(frac_inv(:,9:11),2)];
%% 
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
%% Display both maximum horizontal span and diff betwen ex and in
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 800]);
%Diff between ex and in shwon with histogram and CDF for L23 L4 and L5
subplot(2,1,1);
yyaxis left;ylabel('Counts');
h = histogram(diffL23fr,8);h.BinWidth = 0.01;box off;hold on;h2 = histogram(diffL4fr,8);h2.BinWidth =  0.01;hold on;h3 = histogram(nonzeros(diffL5fr),8);h3.BinWidth =  0.01;
h.EdgeColor = 'k';h.FaceColor = [0.7 0 0];h2.EdgeColor = 'k';h2.FaceColor = [0 1 0];h3.FaceColor = [0 0 1];
yyaxis right;p1=cdfplot(diffL23fr);hold on;p2=cdfplot(diffL4fr);hold on;p3=cdfplot(nonzeros(diffL5fr));grid off; title('');
ylabel('Cumulative');xlabel('Ex-In');p1.Color=[1 0 0];p2.Color=[0 1 0];p3.Color=[0 0 1];legend('L23', 'L4','L5');%legend box off;
subplot(2,1,2);
scatter(ex_spanh,in_spanh,'ko');set(gcf,'color','w');refline(1,0);xlabel('Excitation');ylabel('Inhibition');title('Max. horizontal span');
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,[]);
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
%group based on pial depth
g1=find(pia_input>=220);g2=find(pia_input<220);g3=[];
gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
gv(g1)=1;gv(g2)=2;gv(g3)=3;
%call function
[stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
g1=[];g2=[];g3=[];
% %% Display fraction for ex and in as well as diff for all 16 rows and columns 
% %group based on pial depth
% g1=1:45;g2=46:148;g3=[];
% gv=NaN*ones(1,size(g1,1)+size(g2,1)+size(g3,1));
% gv(g1)=1;gv(g2)=2;gv(g3)=3;
% %call function
% [stats_g] = display_inputs([frac_exv_m frac_inv],[frac_exh frac_inh],frac_diffv,frac_diffh,gv);
% g1=[];g2=[];g3=[];
%% Display correlation between all PCs and pial depth Multiple comparison 
com=[score_ex(:,1:3) score_in(:,1:3) L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23fr diffL4fr pia_input]; 
correlation_matrix(com,1);title('Vertical');
com=[score_ex(:,1:3) score_in(:,1:3) frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' pia_input]; 
correlation_matrix(com,1);title('Horizontal');
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
%IN
corr_plot(L23fr(:,2),pia_input,[],{'L23in input','Pial depth'});set(gca,'Ydir','reverse');
corr_plot(L4fr(:,2),pia_input,[],{'L4in input','Pial depth'});set(gca,'Ydir','reverse');
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
%% Max density per cell apical and basal
max_densba=[max(max(ba_map(:,:,:)))];
max_densba=reshape(max_densba,1,148);
max_densap=max(max(ap_map(:,:,:)));
max_densap=reshape(max_densap,1,148);
%% split in medial and lateral
for i=1:length(nan_vector)
MLba_diff(i)=(sum(nonzeros(ba_map(:,1:8,i)))/length(nonzeros(ba_map(:,1:8,i))))-(sum(nonzeros(ba_map(:,9:end,i)))/length(nonzeros(ba_map(:,9:end,i))));
MLap_diff(i)=(sum(nonzeros(ap_map(:,1:8,i)))/length(nonzeros(ap_map(:,1:8,i))))-(sum(nonzeros(ap_map(:,9:end,i)))/length(nonzeros(ap_map(:,9:end,i))));
end
%% Morphology parameters vs PCs
com=[];com=[score_ex(:,1:3) score_in(:,1:3) morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff']; 
correlation_matrix(com,1);title('Morphology vs PC input scores');
%% %% Morphology parameters vs real parameters vertical
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) diffL23fr diffL4fr morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff']; 
correlation_matrix(com,1);title('Morphology vs Vertical input');
%% %% Morphology parameters vs real parameters horizontal
com=[];com=[frh_medial(:,1)  frh_medial(:,2) frh_lateral(:,1) frh_lateral(:,2) frh_diff_medial frh_diff_lateral (ex_spanh-in_spanh)' morph_parameters(:,1:21) max_densba' max_densap' MLba_diff' MLap_diff']; 
correlation_matrix(com,1);title('Morphology vs Horizontal input');
%% %%Display desired correlations between pial depth and L23/L4
corr_plot(L4fr(:,1),morph_parameters(:,8),morph_parameters(:,9),{'L4ex input','Apical width/height','XSA'});
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
 [od_out sftf_out spon_out] = concat_invivo(L23_PC);
%%  
com=[od_out(:,[1 2 3 4 5 8]) sftf_out spon_out]; 
correlation_matrix(com,1);title('In vivo alone');
 %% %%Display desired correlations between pial depth and L23/L4
corr_plot(com(:,9),com(:,6),com(:,6),{'L4ex input','Apical width/height','XSA'});