%% Clean up
clearvars
close all
% addpath(genpath('C:\Users\drguggiana\Dropbox\Bonhoeffer_code'))
%% Load the gaussian models

%define the main file path
model_path = 'R:\Share\Simon\Drago_Volker_Simon\Trace_mapmaker_out';

%load the file
%define the file tags
file_tags = {'_EImapsSetup1','EImapsSetup2'};
%allocate memory to store the data
all_data = cell(length(file_tags),1);

%for all the tags
for tags = 1:length(file_tags)
    %load the most recent file from the overlap ones
    listing = dir(strcat(model_path,'\*',file_tags{tags},'*.mat'));
    dates = datetime({listing.date});
    [~,ind] = max(dates);
    load_name = listing(ind).name;
    
    %load all the variables in the file
    all_data{tags} = load(fullfile(model_path,load_name));
    all_data{tags} = all_data{tags}.all_maps;
end

%concatenate the data from both setups
all_data = cat(1,all_data{:});
%get the number of cells
cell_num = size(all_data,1);
%% Assemble the data structure to store the entire data set

%create the fields corresponding to the info in the mapmaker files
invitro_struct = cell2struct(all_data,{'cellName','excMap','inhMap','excCounts','inhCount','somaCenter'},2);
%% Load the flip map (cell name, hemisphere(not relevant), slice orientation,
%pial distance, cell ID (from Simon))

%define the map path
excinh_path = 'R:\Share\Simon\Drago_Volker_Simon\';

%load the file
load(fullfile(excinh_path,strcat('Flip_map','.mat')))

%make sure the names of the cells are ordered accordingly in both the file
%and the dataset (output of intersect is sorted by default, so order should be
%compatible with the layers file too)
[~,ia,ib] = intersect({invitro_struct(:).cellName},Flip_map(:,1));
invitro_struct = invitro_struct(ia);
Flip_map = Flip_map(ib,:);
%% Load the layer info for each map

%define the search path
layers_path = 'R:\Share\Simon\Drago_Volker_Simon\layer_GUI_out\';

%generate the possible paths from this main folder
layers_all = strsplit(genpath(layers_path),';')';
layers_all = layers_all(2:3);

%allocate memory to store the matrices
layers_persetup = cell(length(layers_all),1);
%for all the folders
for folders = 1:length(layers_all)
    %get the list and number of files
    file_list = dir(layers_all{folders});
    file_list = file_list(3:end);
    file_num = length(file_list);
    %allocate memory to store the layer maps
    layer_cell = cell(file_num,2);
    
    %load the files (both name and data)
    %for all the files
    for files = 1:file_num
        %get the file name and store it (without extension)
        layer_cell{files,1} = file_list(files).name(1:end-4);
        %and load the data
        layer_cell{files,2} = load(fullfile(file_list(files).folder,file_list(files).name),'new_grid');
        layer_cell{files,2} = layer_cell{files,2}.new_grid;
        
    end
    %store the matrix in the larger cell
    layers_persetup{folders} = layer_cell;
end

%concatenate the two setups
layers_all = cat(1,layers_persetup{:});

%make sure the names of the cells are ordered accordingly in both the file
%and the dataset (output of union is sorted by default, so order should be
%compatible with the layers file too)
[~,ia,ib] = intersect({invitro_struct(:).cellName},layers_all(:,1));
invitro_struct = invitro_struct(ia);
layers_all = layers_all(ib,:);
%% Add the flip and extra info from that file to the structure

%for all the cells
for cells = 1:cell_num
    %add the extra info
   
    invitro_struct(cells).hemisphere = Flip_map{cells,2};
    invitro_struct(cells).sliceOri = Flip_map{cells,3};
    invitro_struct(cells).pialD = Flip_map{cells,4};
    invitro_struct(cells).cellID = Flip_map{cells,5};
   
    invitro_struct(cells).layers = layers_all{cells,2};
end
%sort the cells according to index (so they end up chronologically also)
[~,idx] = sort([invitro_struct(:).cellID]);
invitro_struct = invitro_struct(idx);
%% Apply the flipping of the maps that were oriented with the slide inverted

%run through the maps, flip them if the value is 0
for cells = 1:cell_num
    %if the slice was flipped
    if invitro_struct(cells).sliceOri==0
        %flip the map and layers left to right
        invitro_struct(cells).excMap = invitro_struct(cells).excMap(:,16:-1:1);
        invitro_struct(cells).inhMap = invitro_struct(cells).inhMap(:,16:-1:1);
        invitro_struct(cells).layers = invitro_struct(cells).layers(:,16:-1:1);
    end
end
%% Calculate the excitation and inhibition per layer (normalized and total)

%define the numbers of the target layers
layer_list = [2 3 4 5];
%get the number of layers
layer_num = length(layer_list);
%for all the cells
for cells = 1:cell_num
    %allocate memory to store the fraction vector (including space for exc and inh)
    excinh_perlayer = zeros(layer_num,2);
    %allocate memory for total charge per layer
    excinh_perlayer_total = zeros(layer_num,2);
    %get the layer map
    lyr = invitro_struct(cells).layers;
    %also the exc and inh
    exc = invitro_struct(cells).excMap;
    inh = invitro_struct(cells).inhMap;
    
    %for each layer, calculate the charge normalized to whole cell
    for layers = 1:layer_num
        
        excinh_perlayer(layers,1) = nansum(abs(exc(lyr==layer_list(layers))))/nansum(abs(exc(:)));
        excinh_perlayer(layers,2) = nansum(abs(inh(lyr==layer_list(layers))))/nansum(abs(inh(:)));
        
        excinh_perlayer_total(layers,1) = nansum(abs(exc(lyr==layer_list(layers))));
        excinh_perlayer_total(layers,2) = nansum(abs(inh(lyr==layer_list(layers))));
    end
    %add this information as a field in the structure
    invitro_struct(cells).excFracPerLayer = excinh_perlayer(:,1);
    invitro_struct(cells).inhFracPerLayer = excinh_perlayer(:,2);
    invitro_struct(cells).excRawPerLayer = excinh_perlayer_total(:,1);
    invitro_struct(cells).inhRawPerLayer = excinh_perlayer_total(:,2);
    %also add the total charge per map
    invitro_struct(cells).excinhTotal = [nansum(abs(exc(:)));nansum(abs(inh(:)))];
    
end
%% Calculate side bias per layer, separately for Exc and Inh

close all
%center the cell maps via interpolation based on the soma information

%create a waitbar
w_bar = waitbar(0,'Calculating side bias');
%define the amount of interpolation (the whole grid is 1035umx1035um, with
%69um between each pair of stimulation points)
int_amount = 70;
% %and also the final amplification factor (i.e. size of the map for display)
% map_amount = 20;

%define the field names for the exc and inh map
fields = {'excMap','inhMap'};
%for all the cells
for cells = 1:cell_num
    %show progress
    waitbar(cells/cell_num,w_bar)
    %allocate memory to store the vector (including space for exc and inh)
    sideBias_perlayer = zeros(layer_num,2);
    %get the layer map
    lyr = invitro_struct(cells).layers;
    %get the x soma position of the cell
    soma_x = invitro_struct(cells).somaCenter(1);
    %use imresize for the interpolation and centering of the layer map
    interpol_lyr = round(imresize(lyr,int_amount));
    interpol_lyr = circshift(interpol_lyr,-int16(soma_x).*int_amount/(int_amount-1),2);
    
    %for polarity
    for polars = 1:2
        %get the exc or inh map
        curr_map = invitro_struct(cells).(fields{polars});
        %use imresize for the interpolation
        interpol_map = imresize(curr_map,int_amount);
        %then shift the matrix to put the soma center in the center of the
        %image, and erase the portion that shifted around
        cent_curr = circshift(interpol_map,-int16(soma_x).*int_amount/(int_amount-1),2);
        %turn NaNs into 0
        cent_curr(isnan(cent_curr)) = 0;
        %for every layer
        for layers = 1:layer_num
            %calculate the side bias index per layer
            left_side = double(cent_curr(:,1:size(cent_curr,2)/2)~=0&interpol_lyr(:,1:size(cent_curr,2)/2)==layer_list(layers));
            right_side = double(cent_curr(:,1+size(cent_curr,2)/2:end)~=0&interpol_lyr(:,1:size(cent_curr,2)/2)==layer_list(layers));
            sideBias_perlayer(layers,polars) = (sum(right_side(:)) - sum(left_side(:)))/...
                (sum(right_side(:)) + sum(left_side(:)));
        end
        
    end
    
    %save the indexes in the structure
    invitro_struct(cells).excSideBias = sideBias_perlayer(:,1);
    invitro_struct(cells).inhSideBias = sideBias_perlayer(:,2);
end
%close the waitbar
close(w_bar)
%% Calculate overlap per layer

%for all the cells
for cells = 1:cell_num
    %get the binarized maps
    bin_exc = invitro_struct(cells).excMap;
    bin_exc(isnan(bin_exc)) = 0;
    bin_exc = bin_exc~=0;
    bin_inh = invitro_struct(cells).inhMap;
    bin_inh(isnan(bin_inh)) = 0;
    bin_inh = (bin_inh~=0);
    %get the layers
    lyr = invitro_struct(cells).layers;
    %get the overlap map
    ov_map = squeeze(sum(cat(3,bin_exc,bin_inh),3));
    %allocate a vector to store the overlap per layer
    ov_perlayer = zeros(layer_num,1);
    %for all the layers
    for layers = 1:layer_num

        %get the layers of interest
        ov_layers = ov_map(lyr==layer_list(layers));
        
        %get the index with the corresponding layers
        ov_perlayer(layers) = sum(ov_layers(:)>1)./sum(ov_layers(:)>0);

    end
    %load the info in the final structure
    invitro_struct(cells).overlapPerLayer = ov_perlayer;
end
%% Save the structure

%define the output path
save_path = 'R:\Share\Simon\Drago_Volker_Simon\Full_data_structure';

%assemble the file name
save_name = strcat(datestr(now,'yymmdd_HHMM'),'_dataStruct.mat');

%save the file
save(fullfile(save_path,save_name),'invitro_struct')