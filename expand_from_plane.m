function added_cell = expand_from_plane(str,tc_threshold,selectivity_threshold,rotation_offset,varargin)

%% Get the cell data

% get the number of cells
cell_num = size(str,1);
skip_count = 1;
selected_count = 1;

 % get the fit_resp
 fit_resp = {str.fit_resp};
 fit_vector = zeros(size(fit_resp,1),1);
 for cells = 1:length(str)
     if ~isnan(fit_resp{cells})
        fit_vector(cells) = 1;
     end
 end
 tc_all = cat(2,fit_resp{fit_vector==1});
 
 % get all of the tuning curves
%  tc_all = cat(2,str(~isnan(.fit_resp);

% allocate memory to accumulate the coordinates of the different maps
maps_cell = cell(cell_num,4);
% for all the cells
for cells = 1:cell_num
    if str(cells).iviv == 0 || isempty(str(cells).ORIpref)% || isempty(str(cells).ang_wmap)
        continue
    end
    

    % if the random flag is on, shuffle the angles
    if ~isempty(varargin)
       if varargin{1} == true
          
           % select one at random
           tc = tc_all(:,randperm(size(tc_all,2),1));
       else
           % get the tuning curve
           tc = str(cells).fit_resp;
       end
    else
        % get the tuning curve
        tc = str(cells).fit_resp;
    end
    
    % get the max
    tc_max = max(tc);
    
    if tc_threshold(cells) == 0 
        skip_count = skip_count + 1;
        continue
    else
%         fprintf(strjoin({'Selected cells',num2str(selected_count),'\r\n'},'_'))
        selected_count = selected_count + 1;
    end
    
    % get the excitation and inhibition map
    excMap = str(cells).excMap;
    excMap = excMap./min(excMap(:));
    inhMap = str(cells).inhMap;
    inhMap = inhMap./max(inhMap(:));
    
    
    apiMap = str(cells).morphoMap_apical_aligned;
    basMap = str(cells).morphoMap_basal_aligned;
    
    % get the preferred orientation
    rotation_exc = deg2rad(tc_max-rotation_offset);
    rotation_inh = deg2rad(tc_max-rotation_offset);

    if selectivity_threshold(cells) == 0
%         fprintf(strcat('Cells skipped:',num2str(skip_count) ,'\r\n'))
        skip_count = skip_count + 1;
        continue
    end
    
    [x,y,z,map] = map_3dplotter(excMap,-rotation_exc,1,[],1,0);
    maps_cell{cells,1} = [x,y,z,map];
    [x,y,z,map] = map_3dplotter(inhMap,-rotation_inh,1,[],2,0);
    maps_cell{cells,2} = [x,y,z,map];

    % don't plot structure for cells without structure
    if sum(str(cells).morphoMap_apical_aligned) == 0
        continue
    end
    [x,y,z,map] = map_3dplotter(apiMap,-rotation_inh,0.1,[],3,0);
    maps_cell{cells,3} = [x,y,z,map];
    [x,y,z,map] = map_3dplotter(basMap,-rotation_inh,0.1,[],4,0);
    maps_cell{cells,4} = [x,y,z,map];
end

%% Add the maps

% allocate memory for the output
added_cell = cell(4,1);
% for all the map types
for types = 1:4
    % concatenate the data
    cat_map = cat(1,maps_cell{:,types});
    % if there are no maps here, skip
    if isempty(cat_map)
        added_cell{types} = [];
        continue
    end
    % get the unique coordinates (need to turn to string since unique
    % doesn't work well with floats)
    coord_string = num2str(cat_map(:,1:3));
    % get the unique members and their positions
    [unique_string,ia,ic]  = unique(coord_string,'rows');
    % allocate memory to save the results
    added_maps = zeros(length(ia),4);
    % for all the unique coordinates
    for coord = 1:length(ia)
        % add the intensities at this coordinate
        added_maps(coord,4) = sum(cat_map(ic==coord,4));
        % write the corresponding coordinates
        added_maps(coord,1:3) = str2num(unique_string(coord,:));
    end
    
    % save the maps in the output cell
    added_cell{types} = added_maps;
end