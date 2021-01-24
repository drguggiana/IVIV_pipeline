function centered_soma = calculate_centered_soma(str)
% calculate the soma center of the setup 2 neurons

cell_num = length(str);
% get the slice ori cloud centers for setup 2
centroid_vector = zeros(cell_num,1);
% get the map centroids
for cells = 1:cell_num
    % get the map
    map = str(cells).excMap(3:5,:);
    % get the components
    cc = bwconncomp(map);
    rp = regionprops(cc,map,{'Area','WeightedCentroid'});
    % if there's more than one, leave the largest
    if cc.NumObjects > 1
        areas = cat(1,rp.Area);
        [~,idx] = max(areas);
        rp = rp(idx);
    end
    % store the centroid
    centroid_vector(cells) = rp.WeightedCentroid(:,1);
end
% get the sliceori
slice_ori = cat(1,str.sliceOri);
% get the setup
setup = (1:cell_num)';
setup = setup>47;
% get the mean centers
centered_soma = [mean(centroid_vector(slice_ori==0&setup==1)),...
    mean(centroid_vector(slice_ori==1&setup==1))];