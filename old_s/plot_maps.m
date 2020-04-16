function plot_maps(str,exmap,inhmap,nan_cells,c_id,pia_input)

%% Plot average maps
F = figure;
set(gcf,'color','w');
set(F, 'Name', 'Correlation pial depth');
set(F, 'Position', [200, 0, 800, 200]);
subplot(1,3,1);
sf=1;    
exc_map =exmap(:,:,nan_cells(c_id));
inh_map =inhmap(:,:,nan_cells(c_id))
ove_map = cat(3,exc_map,inh_map);
%define the plot type (2 for excitatory)
explot_type = 2;
if str(c_id).sliceOri==0
   xco=-str(nan_cells(c_id)).somaCenter(1);
else
     xco=-str(nan_cells(c_id)).somaCenter(1);
end
     % p_i=[xco pia_input(c_id)]; 
     p_i=[xco str(nan_cells(c_id)).somaCenter(2)];
     
map_plot3(exc_map,'',explot_type,F,sf,0,1);
hold on;        
      x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
   
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
    plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
   
    hold on;

  subplot(1,3,2);
  inplot_type = 3;
 map_plot3(inh_map,'',inplot_type,F,sf,0,1);
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
   plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5)
    hold on;

  subplot(1,3,3);
   bplot_type = 1;
map_plot3(ove_map,'',bplot_type,F,sf,0,1)
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   

  
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
    plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5)
    hold on;
end
 

