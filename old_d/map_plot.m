function map_plot(data_in,varargin)

%if the 
figure

if nargin > 2
    switch varargin{2}
        case 0 %standard plot the matrix directly
            imagesc(data_in)
            colorbar
        case 1 %plot the E/I overlap
%             %binarize the maps
%             bin_map = data_in>0;
            bin_map = data_in;
  
            %define the smooth factor
            sf = 1;
            %generate a blank matrix to fill in the other color channels
            blank = ones(size(bin_map,1),size(bin_map,2));
            
            %set the excitation map as a red image, smoothing by sf
            exc_map = imresize(cat(3,blank,1-normr_2(bin_map(:,:,1)),1-normr_2(bin_map(:,:,1))),sf);
            %and the inhibition map as a blue one
            inh_map = imresize(cat(3,1-normr_2(bin_map(:,:,2)),1-normr_2(bin_map(:,:,2)),blank),sf);

            % %blend the two images using alpha (and make it double cause default is 8bit
            im_ex = double(imfuse(exc_map,inh_map,'method','blend'));
           
            
            %normalize the image
            im_ex = normr_2(im_ex);
            %remove the NaNs (if there was an empty channel)
            im_ex(isnan(im_ex)) = 0;
            %plot the image
            imagesc(im_ex)
    end
else 
    %standard plot the matrix directly
    imagesc(data_in)
    colorbar
end

%add the lines and labels for the layers
hold('on')
set(gca,'YTick',[1.5, 4.5, 8, 13],'YTickLabels',{'L1','L2/3','L4','DL'},...
    'TickLength',[0 0],'XTick',[])
plot(0:17,2.5.*ones(1,18),'k-')
plot(0:17,6.5.*ones(1,18),'k-')
plot(0:17,9.5.*ones(1,18),'k-')

if nargin > 1
    title(varargin{1})
end