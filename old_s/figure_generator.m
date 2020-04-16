%Load structure with 147 cells used for the paper
str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% Define array names and get the arrays 
field_list = {'subpixel_excMap','subpixel_inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
iviv_cells = find([str(:).iviv]==1);
morpho_cells = ~cellfun(@isempty, {str.morph});
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
% copy the cell to have the original maps later
original_maps = cell_cell;
original_maps_raw = cell_cell_raw;
%Pia vector for ex and inh maps 
pia_input=original_maps(:,end);
%% %% Read out important parameters
ex_map = reshape(original_maps(:,1:256)',16,16,length(str));
in_map = reshape(original_maps(:,257:512)',16,16,length(str));
% Get 16x16 maps for ex and in RAW
ex_map_raw = reshape(original_maps_raw(:,1:256)',16,16,length(str));
in_map_raw = reshape(original_maps_raw(:,257:512)',16,16,length(str));
%Calculate simple difference between maps
diff_map=ex_map-in_map;
% Morphology and Cell ID 
for i=1:length(str)
    cellID_str(i)=str(i).cellID;
    if ~isnan(str(i).morph)==1;
        morph_cells(i)=1;
        morph_parameters(i,:)=str(i).morph;
    else
        morph_cells(i)=0;
        morph_parameters(i,:)=ones(1,24)*NaN;
    end
end
morph_cells_id=find(morph_cells==1);
% Setup A and Setup B
setups=[zeros(47,1);ones(100,1)];
% Slice orientation
slice_ori=[str(:).sliceOri];
%Fractions
frv=reshape([str(:).frac_vert],32,length(str))';
frh=reshape([str(:).frac_hori],32,length(str))';
L23fr=[nanmean(frv(:,3:5),2) nanmean(frv(:,19:21),2)];
L4fr=[nanmean(frv(:,6:7),2) nanmean(frv(:,22:23),2)];
L5fr=[nanmean(frv(:,8:10),2) nanmean(frv(:,24:26),2)];
%Differences of fractions
L5frt=L5fr;
L5frt(find(L5fr(:,1)==0),1)=NaN ;
L5frt(find(L5fr(:,2)==0),2)=NaN ;
diffL23fr=L23fr(:,1)-L23fr(:,2);
diffL4fr=L4fr(:,1)-L4fr(:,2);
diffL5fr=L5frt(:,1)-L5frt(:,2);
%Principal compponent scores
scores=reshape([str(:).PCs],6,length(str))';
%% Calculate the maximum horizontal span overall per layer
[ex_spanhL23 ex_spanhL4 ex_spanhL5 in_spanhL23 in_spanhL4 in_spanhL5] = span_perLayer(ex_map,in_map,nan_vector)
%% Display fraction for ex and in as well as diff for all 16 rows and columns 
[stats_g] = display_inputs([frv],[frh],frv(:,1:16)-frv(:,17:end),frh(:,1:16)-frh(:,17:end),[]);