function [map_out,curr_info_out] = interpolate_map(curr_info,curr_map)

%get the info for the current map
%     curr_info = trace2folder(trace2folder(:,2)==maps,:);
% %get the polarity
% polarity = mode(curr_info(:,3));

num_positions = size(curr_map,2);
%     %if an inhibitory map, skip it
%     if polarity == 1
%         continue
%     end

% copy the info
curr_info_out = curr_info;

%check if there are any direct responses
num_direct = sum(curr_info(:,5)==0);
%if no direct responses, skip the interpolation and load the binned map
if num_direct == 0
    map_out = curr_map;
    return
end

%     %don't interpolate the post TTX maps
%     if mode(time_pertrace(trace2folder(:,2)==maps))==0
%         continue
%     end

%blank the direct traces
curr_map(:,curr_info(:,5)==0) = NaN;

%reshape the current map to a 3d grid (i.e. x , y and time)
curr_map = reshape(curr_map,[],sqrt(num_positions),sqrt(num_positions));

%set up the interpolation

%get the 3d coordinates of the direct responses
direct_coord = find(isnan(curr_map));
[qx,qy,qz] = ind2sub(size(curr_map),direct_coord);
%and of the non-direct responses
nondirect_coord = find(~isnan(curr_map));
[dx,dy,dz] = ind2sub(size(curr_map),nondirect_coord);


filled_coord = griddatan([dx,dy,dz],curr_map(nondirect_coord),[qx,qy,qz],'linear');

filled_map = curr_map;
filled_map(direct_coord) = filled_coord;
%store the interpolated map
map_out = reshape(filled_map,[],num_positions);
% also correct the info to change the direct traces into not
curr_info_out(curr_info_out(:,5)==0,5) = 4;