function [aligned_map,aligned_soma] = align_subpixel(map,slice_ori,soma_slice,soma_center,setup)
% align a map to a precalculated center based on subpixel interpolation

% define the grid spacing in microns
grid_spacing = 69;

% get the map size in microns
map_size = round(16.*grid_spacing);
% get the map limits
map_lim = map_size/2-grid_spacing/2;

% create the grid
[Y,X] = ndgrid(-map_lim:grid_spacing:map_lim,...
    -map_lim:grid_spacing:map_lim);


% select the appropriate offset
if slice_ori == 0
    offset = (soma_slice(1)-8.5)*grid_spacing;
    % add the corrected soma position
    aligned_soma = [-soma_center(1),...
        soma_center(2)];
else
    offset = (soma_slice(2)-8.5)*grid_spacing;
    % add the corrected soma position
    aligned_soma = [soma_center(1),...
        soma_center(2)];
end
% subtract the offset
center_X = X + offset;

% if it's a setup 1 cell, don't interpolate it
if setup == 0
    aligned_map = map;
else
    % allocate memory for the output
    aligned_map = zeros(size(map));
    % for all the maps in map
    for maps = 1:size(map,3)
        % interpolate the map
        interpolant = griddedInterpolant(Y,X,map(:,:,maps));
        % replace the map
        aligned_map(:,:,maps) = interpolant(Y,center_X);
    end

%     % save the raw map also for examples
%     aligned_map = new_map;
end
