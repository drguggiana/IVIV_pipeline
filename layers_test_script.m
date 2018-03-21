close all
%define the image dimensions in microns
image_xmicro = 2548;
image_ymicro = 1903;
%define the number of points in the grid
x_points = 16;
y_points = 16;
%load the image
% im_load = imread('I:\Simon Weiler\EXPLORER ONE\161013iviv\SW0002\images\videoImg_SW0002_1.tif'); 
im_load = imread('I:\Simon Weiler\EXPLORER ONE\180219iviv\i1.tif');
%get the image dimensions in pixels
[image_ypix,image_xpix] = size(im_load);

%calculate the micron to pixel conversion factor
micro2pix_x = image_xpix/image_xmicro;
micro2pix_y = image_ypix/image_ymicro;

%plot it
imagesc(im_load)
colormap(gray)
axis equal
% %send the image to appdata
% setappdata(gcf,'im_load',im_load)

%load the grid data
% grid_path = 'I:\Simon Weiler\EXPLORER ONE\161013iviv\SW0002\map01\SW0002MAAA0001.xsg';
grid_path = 'I:\Simon Weiler\EXPLORER ONE\180219iviv\SW0002\map01\SW0002MAAA0001.xsg';


xsg_data = load(grid_path,'-mat');
%load the grid spacing (convert to pixels)
grid_xspace = xsg_data.header.mapper.mapper.xSpacing*micro2pix_x;
grid_yspace = xsg_data.header.mapper.mapper.ySpacing*micro2pix_y;

%load the grid offset (convert to pixels)
grid_xoffset = xsg_data.header.mapper.mapper.xPatternOffset*micro2pix_x;
grid_yoffset = xsg_data.header.mapper.mapper.yPatternOffset*micro2pix_y;
% grid_xoffset = xsg_data.header.mapper.mapper.xPatternOffset;
% grid_yoffset = xsg_data.header.mapper.mapper.yPatternOffset;

%load the grid rotation (convert to pixels)
grid_rotation = xsg_data.header.mapper.mapper.spatialRotation;

%calculate the coordinates of the grid vertices (0 rotation and offset
%first)
%allocate memory for the coordinates
grid_coord = zeros(y_points*x_points,2);
%also for a matrix map of the points
grid_map = zeros(y_points,x_points);
%start a counter for overall positions
pos_count = 1;
%for all the x positions
for x = 1:x_points
    %for all the y positions
    for y = 1:y_points
        %fill up the coordinates based on the respective spacings
        grid_coord(pos_count,:) = [1+grid_yspace*(y-1),1+grid_xspace*(x-1)];
        %also store the index of the point in the corresponding map
        %position
        grid_map(y,x) = pos_count;
        %update the counter 
        pos_count = pos_count + 1;
    end
end
%center the grid around 0 to rotate
%calculate the grid original center
y_center = (max(grid_coord(:,1)) - min(grid_coord(:,1)))/2 + 1;
x_center = (max(grid_coord(:,2)) - min(grid_coord(:,2)))/2 + 1;

grid_coord(:,1) = grid_coord(:,1) - y_center;
grid_coord(:,2) = grid_coord(:,2) - x_center;

% Create rotation matrix
R = [cosd(grid_rotation) -sind(grid_rotation); sind(grid_rotation) cosd(grid_rotation)];
%rotate the array
grid_coord = (R*grid_coord')';
%return it to its original position
grid_coord(:,1) = grid_coord(:,1) + y_center;
grid_coord(:,2) = grid_coord(:,2) + x_center;
%grid is now at the top corner shift it to the image center
%calculate the distances from the edge of the grid to the end (accounting
%for the centering above)
% y_toend = image_ypix-max(grid_coord(:,1));
% x_toend = image_xpix-max(grid_coord(:,2));

y_toend = image_ypix/2 - y_center;
x_toend = image_xpix/2 - x_center;

%shift by half that distance
grid_coord(:,1) = grid_coord(:,1) + y_toend - grid_yoffset;
grid_coord(:,2) = grid_coord(:,2) + x_toend + grid_xoffset;

%plot the grid on the image
hold('on')
plot(grid_coord(:,2),grid_coord(:,1),'*')
plot(image_xpix/2,image_ypix/2,'ro')
% plot(x_center,y_center,'ro')
% set(gca,'XLim',[-100 100],'YLim',[-100 100])

% %apply offset and rotation
% grid_coord2 = grid_coord;
% grid_coord2(:,1) = grid_coord2(:,1) + ;
% grid_coord2(:,2) = grid_coord2(:,2) + ;