function map_plot3(data_in,varargin)

%arguments
%1) data
%2) label (default no label)
%3) plot matrix directly (0, default) vs the overlap plot (1)
%4) target figure handle (default none, create new figure)
%5) smoothing factor (default 1)
%6) if 1, draw crosshair in map

%if the fourth argument is specified, plot in the provided figure handle
if nargin > 3
    figure(varargin{3})
else %otherwise create a new figure
    figure
end

%if a smoothing factor was specified
if nargin > 4
    sf = varargin{4};
    
else% otherwise, assume 1
    sf = 1;
end
%also define the scaling factor for the plot labels, combining the input
%smoothing factor plus the size of the image
sf_plot = sf*size(data_in,1)/16;

%if the third argument (plot type) is specified
if nargin > 2
    switch varargin{2}
        case 0 %standard plot the matrix directly
            imagesc(data_in)
            colorbar
        case 1 %plot the E/I overlap


            bin_map = abs(data_in)>0;

  
            %generate a blank matrix to fill in the other color channels
            blank = ones(size(bin_map,1),size(bin_map,2));
            
            %set the excitation map as a red image, smoothing by sf
            exc_map = imresize(cat(3,blank,1-normr_2(bin_map(:,:,1)),1-normr_2(bin_map(:,:,1))),sf);

            %and the inhibition map as a blue one
            inh_map = imresize(cat(3,1-normr_2(bin_map(:,:,2)),1-normr_2(bin_map(:,:,2)),blank),sf);

            
            % %blend the two images using alpha (and make it double cause default is 8bit
            im_ex = double(imfuse(exc_map,inh_map,'method','blend'));
            
            %restore the NaNs
            nan_map = imresize(squeeze(sum(isnan(data_in),3)),sf);
            im_ex(cat(3,nan_map,nan_map,nan_map)>0) = NaN;
           
           
            %normalize the image
            im_ex = normr_2(im_ex);

            %set the NaNs alpha to 0 so they show the background
            imAlpha = ones(size(im_ex,1),size(im_ex,2));
            imAlpha(isnan(sum(im_ex,3))) = 0;
            
            %plot the image
            image(im_ex,'AlphaData',imAlpha)
            
            color_column = ((0:255)/255)';
            color_column2 = ((255:-1:0)/255)';
            colormap([color_column,zeros(256,1),color_column2])
            colorbar(gca,'Ticks',[0 1],'TickLabels',{'Inh','Exc'})
            set(gca,'Color',[0 0 0])

    end
else 
    %standard plot the matrix directly
    imagesc(data_in)
    colormap
end

%add the lines and labels for the layers
%(mode layer assignment: L1:1,2 L2/3:3,4,5,6 L4:7,8 L5:9,10,11 L6:12,13,14
%WM:15,16 )
hold('on')
set(gca,'YTick',[1.5, 4.5, 7.5, 10, 13, 15.5].*sf_plot,'YTickLabels',{'L1','L2/3','L4','L5','L6','WM'},...
    'TickLength',[0 0],'XTick',[])
plot(linspace(0,17*sf_plot,18),2.5.*ones(1,18).*sf_plot,'k-')
plot(linspace(0,17*sf_plot,18),6.5.*ones(1,18).*sf_plot,'k-')
plot(linspace(0,17*sf_plot,18),8.5.*ones(1,18).*sf_plot,'k-')
plot(linspace(0,17*sf_plot,18),11.5.*ones(1,18).*sf_plot,'k-')
plot(linspace(0,17*sf_plot,18),14.5.*ones(1,18).*sf_plot,'k-')
% plot(0:17,6.5.*ones(1,18),'k-')
% plot(0:17,9.5.*ones(1,18),'k-')

%if a 6th argument is provided, plot a cross in the center of the map
if nargin > 5

    %get the plot limits
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    %get the plot centers
    x_cent = sum(x_lim)/2;
    y_cent = sum(y_lim)/2;
    plot([x_cent x_cent],y_lim,'g-')
    plot(x_lim,[y_cent y_cent],'g-')
end

if nargin > 1
    title(varargin{1})
end