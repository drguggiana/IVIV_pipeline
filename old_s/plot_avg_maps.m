%% Plot average maps for all cells
function plot_avg_maps(str,nan_vector,ex_map,in_map,pia_input,sf,clu,idx_input);
if clu==0
F = figure;
set(gcf,'color','w');
set(F, 'Name', 'Correlation pial depth');
set(F, 'Position', [200, 0, 800, 200]);
subplot(1,3,1);
sf=sf;    
exc_map =nanmean(ex_map(:,:,nan_vector),3);
inh_map =nanmean(in_map(:,:,nan_vector),3);
ove_map = cat(3,exc_map,inh_map);
%define the plot type (2 for excitatory)
explot_type = 2;
%get the pial distance
soma_info = pia_input(1);
map_plot3(exc_map,'',explot_type,F,sf,0,1);
hold on;        
      x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
    %plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3);
   
    hold on;
  end
  subplot(1,3,2);
  inplot_type = 3;
 map_plot3(inh_map,'',inplot_type,F,sf,0,1);
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
   % plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3)
    hold on;
  end
  subplot(1,3,3);
   bplot_type = 1;
map_plot3(ove_map,'',bplot_type,F,sf,0,1)
  x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');   
  for i=1:length(nan_vector)
      xco=str(nan_vector(i)).somaCenter(1)
      p_i=[xco pia_input(i)]; 
    adj_x = ((16*69/2)-p_i(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y = ((16*69/2)-p_i(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
   % plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3)
    hold on;
  end
  
  
  
else
    %Clusters
  %% Ex and in avergae maps
F = figure;
set(F, 'Name', 'Heat average input maps for cluster');
set(F, 'Position', [0, 0, 1000, 1000]);set(gcf,'color','w')
 set(gcf,'color','w');
sf=sf;
clu_num=max(idx_input);
for clu = 1:clu_num
    p_i=[0 nanmean(pia_input(find(idx_input==clu)))];   
       excc_map =nanmean(ex_map(:,:,idx_input==clu),3)
       inhc_map =nanmean(in_map(:,:,idx_input==clu),3)
       ove_map = cat(3,excc_map,inhc_map);
         %define the plot type (2 for excitatory)
        explot_type = 2;
         %get the pial distance
        soma_info = pia_input(idx_input==clu);
        subplot(3,4,clu)
        map_plot3(excc_map,'',explot_type,F,sf,1,1);
             title(num2str(sum(idx_input==clu)));
        inplot_type = 3;
        subplot(3,4,clu+4);
        map_plot3(inhc_map,'',inplot_type,F,sf,1,1);
  
   bplot_type = 1;
   subplot(3,4,clu+8);
   map_plot3(ove_map,'',bplot_type,F,sf,1,1);
    p_i=[];
    
end  
end
end