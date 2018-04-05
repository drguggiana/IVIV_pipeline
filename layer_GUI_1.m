function varargout = layer_GUI_1(varargin)
% LAYER_GUI_1 MATLAB code for layer_GUI_1.fig
%      LAYER_GUI_1, by itself, creates a new LAYER_GUI_1 or raises the existing
%      singleton*.
%
%      H = LAYER_GUI_1 returns the handle to a new LAYER_GUI_1 or the handle to
%      the existing singleton*.
%
%      LAYER_GUI_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAYER_GUI_1.M with the given input arguments.
%
%      LAYER_GUI_1('Property','Value',...) creates a new LAYER_GUI_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before layer_GUI_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to layer_GUI_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help layer_GUI_1

% Last Modified by GUIDE v2.5 05-Mar-2018 09:16:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @layer_GUI_1_OpeningFcn, ...
                   'gui_OutputFcn',  @layer_GUI_1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before layer_GUI_1 is made visible.
function layer_GUI_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to layer_GUI_1 (see VARARGIN)

% Choose default command line output for layer_GUI_1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes layer_GUI_1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%put the cell paths in appdata
setappdata(gcf,'im_path',varargin{2})
setappdata(gcf,'xsg_path',varargin{3})
%load the text counters
set(handles.curr_cell,'String',num2str(varargin{4}))
set(handles.tot_cells,'String',num2str(varargin{5}))
%and the image size in microns
setappdata(gcf,'x_micro',varargin{6}(1))
setappdata(gcf,'y_micro',varargin{6}(2))

%load the image
im_load = imread(getappdata(gcf,'im_path'));

%send the image to appdata (only one channel in case it's a bmp with 3)
setappdata(gcf,'im_load',im_load(:,:,1))

%load the xsg data
xsg_data = load(getappdata(gcf,'xsg_path'),'-mat');

%send the xsg data to appdata
setappdata(gcf,'xsg_data',xsg_data)

%configure the progress bar
prog_bar(handles,0,4)
%set the current axes to axis 1
axes(handles.axes1)
%plot the image and grid
draw_image(handles,0)


% --- Outputs from this function are returned to the command line.
function varargout = layer_GUI_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in draw_1.
function draw_1_Callback(hObject, eventdata, handles)
% hObject    handle to draw_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear the axis
cla
%define the number of layers to use
layer_num = 4;

%draw the image and grid
draw_image(handles,0)
%create a position constraint function to not draw outside the window
fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),...
   get(gca,'YLim'));

%allocate memory to store the coordinates
layer_data = cell(layer_num,1);
%for all the layer boundaries
for layers = 1:layer_num
    
    %move the layer marker to the desired position
    prog_bar(handles,layers,layer_num)
    %start the drawing interface
    line_h = impoly(gca,'Closed',false,'PositionConstraintFcn',fcn);

    %get the position of the line
    % line_coord = getPosition(line_h)
    line_coord = wait(line_h);
    
    %store the coordinates in the target cell
    layer_data{layers} = line_coord;
    

    %plot the line
    hold('on')
    plot(line_coord(:,1),line_coord(:,2),'o-')

end
%complete the progress bar
prog_bar(handles,layers+1,layer_num)
%get the grid assignments based on the lines
new_grid = layer_assign(handles,layer_data);
%send the output to appdata
setappdata(gcf,'new_grid',new_grid)
%also send the layers to appdata
setappdata(gcf,'layer_data',layer_data)
%redraw the image with the new grid
cla
draw_image(handles,1)


% --- Executes on button press in done_1.
function done_1_Callback(hObject, eventdata, handles)
% hObject    handle to done_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%export the coordinate data to the main workspace
% layer_data = getappdata(gcf,'layer_data');
% assignin('base','layer_data',layer_data)

new_grid = getappdata(gcf,'new_grid');
assignin('base','new_grid',new_grid)

%close the window
delete(gcf)

function draw_image(handles,flag)

cla
%load the image from appdata
im_load = getappdata(gcf,'im_load');

%load the xsg data
xsg_data = getappdata(gcf,'xsg_data');

%load the image dimensions in microns
image_xmicro = getappdata(gcf,'x_micro');
image_ymicro = getappdata(gcf,'y_micro');
%define the number of points in the grid
x_points = 16;
y_points = 16;
%load the image

%get the image dimensions in pixels
[image_ypix,image_xpix] = size(im_load);

%calculate the micron to pixel conversion factor
micro2pix_x = image_xpix/image_xmicro;
micro2pix_y = image_ypix/image_ymicro;

%plot it
imagesc(im_load)
colormap(gray)
axis equal

%load the grid spacing (convert to pixels)
grid_xspace = xsg_data.header.mapper.mapper.xSpacing*micro2pix_x;
grid_yspace = xsg_data.header.mapper.mapper.ySpacing*micro2pix_y;

%load the grid offset (convert to pixels)
grid_xoffset = xsg_data.header.mapper.mapper.xPatternOffset*micro2pix_x;
grid_yoffset = xsg_data.header.mapper.mapper.yPatternOffset*micro2pix_y;

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
y_toend = image_ypix/2 - y_center;
x_toend = image_xpix/2 - x_center;

%shift by half that distance
grid_coord(:,1) = grid_coord(:,1) + y_toend - grid_yoffset;
grid_coord(:,2) = grid_coord(:,2) + x_toend + grid_xoffset;

%plot the grid on the image
hold('on')
if flag == 0
    %just plot the original grid and the traces on top
    plot(grid_coord(:,2),grid_coord(:,1),'b*')
    
else
    %plot the grid with layer color assignments
    %load the new grid
    new_grid = getappdata(gcf,'new_grid');
    
    %load the grid map
    grid_map = getappdata(gcf,'grid_map');
    %load the layer data
    layer_data = getappdata(gcf,'layer_data');
    
    %get the number of layers
    layer_num = size(layer_data,1);
    %replace the 0s with layer_num + 1 so the non-assigned points show up in white
    new_grid(new_grid==0) = layer_num + 1;
    %create a color map for the layers
    c_map = [jet(layer_num+1);[1 1 1]];
    
    %for all the grid points
    for points = 1:numel(new_grid)
        %plot each point with it's corresponding color
        plot(grid_coord(points,2),grid_coord(points,1),'b*','Color',c_map(new_grid(grid_map==points),:))
    end
    
    %also plot the lines
    %for all the layers
    for layers = 1:layer_num
        plot(layer_data{layers}(:,1),layer_data{layers}(:,2),'o-','Color',c_map(layers,:))
    end
end

%calculate the new grid center and plot
new_x_center = x_center + x_toend + grid_xoffset;
new_y_center = y_center + y_toend - grid_yoffset;
plot(new_x_center,new_y_center,'ro')
%get the soma coordinates
soma_coord = xsg_data.header.mapper.mapper.soma1Coordinates .* [micro2pix_y,micro2pix_x];
%coordinates are given with respect to image center, so correct
somax = soma_coord(1) + image_xpix/2;
somay = -soma_coord(2) + image_ypix/2;
%plot the center
plot(somax,somay,'go')

%define the spacing in microns of the ticks
tick_space = 200;
%calculate the axes with respect to the grid center
x_axis = (0:image_xmicro) - new_x_center/micro2pix_x;
y_axis = (0:image_ymicro) - new_y_center/micro2pix_y;
%calculate the tick positions
x_ticks = x_axis(1:tick_space:end).*micro2pix_x + new_x_center;
y_ticks = y_axis(1:tick_space:end).*micro2pix_y + new_y_center;
%and define the labels
x_labels = x_axis(1:tick_space:end);
y_labels = -y_axis(1:tick_space:end);

%add axis in microns with respect to the grid center
set(gca,'XTick',x_ticks,'XTickLabels',x_labels,...
    'YTick',y_ticks,'YTickLabels',y_labels,...
    'XTickLabelRotation',45)
xlabel('Distance (um)')
ylabel('Distance (um)')
%send the array points to appdata
setappdata(gcf,'grid_coord',grid_coord)
%and the map
setappdata(gcf,'grid_map',grid_map)
%and plot the traces
plot_traces(handles)


function [new_grid] = layer_assign(handles,layer_data)
%load the array points from appdata
grid_coord = getappdata(gcf,'grid_coord');
%and the map
grid_map = getappdata(gcf,'grid_map');
%create a working copy of the data
grid_copy = grid_coord;
%allocate memory to store the labeled grid
new_grid = zeros(size(grid_map));
%get the number of layers
layer_num = size(layer_data,1);

%for all the layers
for layers = 1:layer_num
    %interpolate the line at all the x positions of the grid points
    line_points = interp1(layer_data{layers}(:,1),layer_data{layers}(:,2),grid_copy(:,2));
    %determine which points are above the line and mark them
    target_points = find(grid_copy(:,1)<line_points);
    %make the points NaN in the coordinates
    grid_copy(grid_copy(:,1)<line_points,2) = NaN;
    %for all the target points
    for points = 1:length(target_points)
        %fill in the assignment matrix
        new_grid(grid_map==target_points(points)) = layers;
    end
    %if it's the last layer, also assign the points below that line (cause
    %actual last layer, not last boundary)
    if layers == layer_num
        %determine which points are above the line and mark them
        target_points = find(grid_copy(:,1)>line_points);
        %for all the target points
        for points = 1:length(target_points)
            %fill in the assignment matrix
            new_grid(grid_map==target_points(points)) = layers+1;
        end
        
    end
end

function [trace_mat] = trace_parser(handles)

%get the trace data from appdata
xsg_data = getappdata(gcf,'xsg_data');

%load the map to get the order of the traces
map_mat = xsg_data.header.mapper.mapper.mapPatternArray;
%get the dimensions of the stimulation array
array_dim = size(map_mat);

%Load the data from the structure
trace_data = xsg_data.data.ephys.trace_1;

%reshape as a 3d matrix
trace_data = reshape(trace_data,[],array_dim(1)*array_dim(2));
% trace_data = permute(trace_data,[2 3 1]);

%define the interval to grab. Each recording is a second and the desired
%range is from 100 to 300 ms
trace_range = 1000:3000;
%trim traces to only the relevant response
trace_data = trace_data(trace_range,:);
%allocate memory to store the reordered data
trace_reorder = zeros([array_dim size(trace_data,1)]);
%and reorder them according to the stimulation randomization
%for all the cells in x
for x = 1:array_dim(1)
    %for all the cells in y
    for y = 1:array_dim(2)
        %get the correct x and y for the current trace
        trace_reorder(x,y,:) = trace_data(:,map_mat(x,y));
    end
end

%output the matrix
trace_mat = trace_reorder;

function plot_traces(handles)

%define the subsampling ratio for time
sub_rat = 80;
%define the amplitude factor
amp_rat = handles.amp_rat_slider.Value;

%load the grid coordinates
grid_coord = getappdata(gcf,'grid_coord');

%process the trace data to plot on top of the pic
trace_mat = trace_parser(handles);
%%
%for all the traces
for x = 1:size(trace_mat,1)
    for y = 1:size(trace_mat,2)
        %get the index corresponding to this position in single index
        curr_ind = sub2ind([size(trace_mat,1),size(trace_mat,2)],x,y);
        %get the x vector, correcting for the array position
        x_vec = (1:length(squeeze(trace_mat(x,y,1:sub_rat:end)))) + grid_coord(curr_ind,2);
        %and the y vector
        y_vec = squeeze(trace_mat(x,y,1:sub_rat:end))./-amp_rat + grid_coord(curr_ind,1);
        %plot the result
        plot(x_vec,y_vec,'r')
    end
end


% --- Executes on slider movement.
function amp_rat_slider_Callback(hObject, eventdata, handles)
% hObject    handle to amp_rat_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
draw_image(handles,0)

% --- Executes during object creation, after setting all properties.
function amp_rat_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amp_rat_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function curr_cell_Callback(hObject, eventdata, handles)
% hObject    handle to curr_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curr_cell as text
%        str2double(get(hObject,'String')) returns contents of curr_cell as a double


% --- Executes during object creation, after setting all properties.
function curr_cell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curr_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tot_cells_Callback(hObject, eventdata, handles)
% hObject    handle to tot_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tot_cells as text
%        str2double(get(hObject,'String')) returns contents of tot_cells as a double


% --- Executes during object creation, after setting all properties.
function tot_cells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tot_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function prog_bar(handles,curr_level,tot_level)

%create the progress bar vector
prog_vec = zeros(tot_level+1,1);
%set the second axis as target
axes(handles.axes2)
%modify the progress bar
prog_vec(1:curr_level) = 1;
imagesc(prog_vec)
%configure the progress bar axes
prog_labels = {'L6','L5','L4','L2/3','L1'};
set(handles.axes2,'YTick',1:5,'YTickLabels',prog_labels)
colormap(handles.axes2,[0.94 0.94 0.94;1 0 0])
%return the target to the first axis
axes(handles.axes1)
