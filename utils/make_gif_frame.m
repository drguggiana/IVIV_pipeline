function make_gif_frame(frame_handle,frame_idx,filename,varargin)
% make a gif

if length(varargin) >= 1
    frametime = varargin{1};
else
    frametime = 0.5;
end

drawnow
% Capture the plot as an image
frame = getframe(frame_handle);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if frame_idx == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',frametime);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',frametime);
end
