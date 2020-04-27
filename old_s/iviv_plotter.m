function iviv_plotter(str,i)
oris=[0:45:315];
%--------------MORPHOLOGY

     ha = tightPlots(1, 4, 6, [1 1], [0.2 0.2], [0.1 0.4], [0.1 0.1], 'centimeters');set(gcf,'color','w');
%   f = [1 5 10 15]; x = 0:0.05:10;
%   for i = 1:length(f)
%       y = sin(f(i) * x);
%       axes(ha(i)); plot(x, y)
%   end
% F=figure;
% set(gcf, 'Position', [300, 300, 300, 125]);
% set(gcf,'color','w');
%H=subplot(1,4,1);

axes(ha(1));
   if ~isnan(str(i).morph(1))==1
       
       m=plot_tree(str(i).morphtraces{1,1},[1 0 0],[0 str(i).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = 'k';
       m1=plot_tree(str(i).morphtraces{1,2},[0 0 0],[0 str(i).morphtraces{1,5} 0],[],1,'-b');hold on;
       m1.EdgeColor = [0.5 0.5 0.5];
       %m1.EdgeColor = 'm';
       m2=plot_tree(str(i).morphtraces{1,3},[0 0 1],[0 str(i).morphtraces{1,5} 0],[],1,'-b');hold on;
       m2.EdgeColor = [0.5 0.5 0.5];
       set(gca,'Ydir','reverse');
        xlim([-255 255]);
 ylim([0 500]);
 axis off;
 set(gcf,'color','w');
 hold on;text(-240,-115,['#' num2str(i)]);
       
else
    m=plot(str(i).pialD,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4);
 set(gca,'Ydir','reverse');
        xlim([-255 255]);
 ylim([0 500]);
 axis off;
 set(gcf,'color','w');
 hold on;text(-240,-115,['#' num2str(i)]);
   end


 %set(H, 'Position', [0.1, 0.3, 0.1, 1]);
  
    
%-----------Polarplot_ipsi_contra
axes(ha(2));
if str(i).iviv==1
if str(i).resp==1
%set(H2, 'Position', [0.2, 0.4, 0.1, 0.3]);

% ax.ThetaGrid = 'off';
% Ax.RGrid = 'off';

    if str(i).contra==1
        polarplot(deg2rad(oris([1:end 1])),str(i).TCpref([1:8 1])/max(str(i).TCpref),'Color',[0 0.4 0.6]);hold on;
        ax = gca;rticks([0:0.5:1]);ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
        ax.RTickLabel = []; 
        ax.ThetaTickLabel = [];
 hold on;title([num2str(round(str(i).OSIpref,2)) ' / ' num2str(round(str(i).DSIpref,2))],'Color',[0 0.4 0.6],'FontWeight','normal');
    elseif str(i).ipsi==1
            polarplot(deg2rad(oris([1:end 1])),str(i).TCpref([1:8 1])/max(str(i).TCpref),'Color',[0.7 0 0]);hold on;
        ax = gca;rticks([0:0.5:1]);ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
        ax.RTickLabel = []; 
        ax.ThetaTickLabel = [];
        hold on;title([num2str(round(str(i).OSIpref,2)) ' / ' num2str(round(str(i).DSIpref,2))],'Color',[0.7 0 0],'FontWeight','normal');
    
    
    else str(i).bino==1
          polarplot(deg2rad(oris([1:end 1])),str(i).TuningCurve([1:8 1])/max(str(i).TuningCurve),'Color',[0 0.4 0.6]);hold on;
polarplot(deg2rad(oris([1:end 1])),str(i).TuningCurve([9:end 9])/max(str(i).TuningCurve),'Color',[0.7 0 0]);
ax = gca;rticks([0:0.5:1]);ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
ax.RTickLabel = []; 
ax.ThetaTickLabel = []; 
         hold on;title([num2str(round(str(i).OSIpref,2)) ' / ' num2str(round(str(i).DSIpref,2))],'Color',[0.5 0.5 0.5],'FontWeight','normal');
    end

 
else
      polarplot(deg2rad(oris([1:end 1])),ones(1,9)*NaN);hold on;
polarplot(deg2rad(oris([1:end 1])),ones(1,9)*NaN);
ax = gca;rticks([0:0.5:1]);ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
ax.RTickLabel = []; 
ax.ThetaTickLabel = [];  
  hold on;title(['nvr'],'FontWeight','normal'); 
    
end
else
  polarplot(deg2rad(oris([1:end 1])),str(i).TuningCurve([1:8 1])/max(str(i).TuningCurve));hold on;
polarplot(deg2rad(oris([1:end 1])),str(i).TuningCurve([9:end 9])/max(str(i).TuningCurve));
ax = gca;rticks([0:0.5:1]);ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
ax.RTickLabel = []; 
ax.ThetaTickLabel = [];  
hold on;title(['n.c.'],'FontWeight','normal');
axis off;
end




 hold on; 
 axes(ha(3));
sf=1
    %plot the excitatory map
    %get the map
    exc_map = str(i).subpixel_excMap;
    %check if the map exists, otherwise skip
    if ~isnan(sum(exc_map(:)))
        %define the plot type (2 for excitatory)
        plot_type = 2;
        %get the pial distance
        soma_info = [str(i).subpixel_soma(1) str(i).subpixel_soma(2)];
        %plot the map
     % H5=subplot(1,4,3);
      %set(H5, 'Position', [0.55, 0.4, 0.2, 0.4]);
        map_plot3_SW(exc_map,'',plot_type, gcf,sf,0,1,soma_info);
         %axis off;
       ax = gca;
       ax.XTickLabel=[];
       ax.YTickLabel=[];
    end
%----------------------------Input map Inhibition------------------------
%subplot(1,7,6);
 hold on;  
 axes(ha(4));
sf=1
    %plot the excitatory map
    %get the map
    inh_map = str(i).subpixel_inhMap;
    %check if the map exists, otherwise skip
    if ~isnan(sum(inh_map(:)))
        %define the plot type (2 for excitatory)
        plot_type = 3;
        %get the pial distance
        soma_info = [str(i).subpixel_soma(1) str(i).pialD];
        %plot the map
      %H6=subplot(1,4,4);
      %set(H6, 'Position', [0.75, 0.4, 0.2, 0.4]);
        map_plot3_SW(inh_map,'',plot_type,gcf,sf,0,1,soma_info);
           ax = gca;
       ax.XTickLabel=[];
       ax.YTickLabel=[];
    end


end