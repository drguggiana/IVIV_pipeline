%% This script plots the figures, for Weiler et al. 2020; USE EVALUATE SECTIONS to go stepy by step through the script
%You will need: data structure str, embedding for UMAP, uipickfiles,
%display_inputs,extent_diff, centroid_plot, corr_plot, correlation_matrix, autoArrangeFigures;
%plotting_embedding_str,
%% START
%Load structure with 147 cells used for the paper
str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str\';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% %% Read out important parameters
% get a vector with the cells to use
iviv_cells = find([str(:).iviv]==1);
%Get normalized maps + raw maps
for i=1:length(str)
ex_map(:,:,i) = str(i).subpixel_excMap;
in_map(:,:,i) = str(i).subpixel_inhMap;
% ex_map_raw(:,:,i) = str(i).subpixel_raw_excMap;
% in_map_raw(:,:,i) = str(i).subpixel_raw_inhMap;
end
%Calculate simple difference between maps
diff_map=ex_map-in_map;
% Morphology and Cell ID 
for i=1:length(str)
    cellID_str(i)=str(i).cellID;
    if ~isnan(str(i).morph)==1;
        morph_cells(i)=1;
        morph_parameters(i,:)=str(i).morph;
    else
        morph_cells(i)=0;
        morph_parameters(i,:)=ones(1,24)*NaN;
    end
end
morph_cells_id=find(morph_cells==1);
%Get the morphology density maps
for i=1:length(str)   
        morpho_basal{i,:}=str(i).morphoMap_basal_aligned;
        morpho_apical{i,:}=str(i).morphoMap_apical_aligned;     
end
for i=1:length(str)
     if isempty(morpho_basal{i,:,:})==0
     ba_map(:,:,i)=morpho_basal{i,:,:};
     ap_map(:,:,i)=morpho_apical{i,:,:};
     else
     ba_map(:,:,i)=ones(16,16)*NaN;
     ap_map(:,:,i)=ones(16,16)*NaN;
     end
end
%Sum for each morpho density map
for i=1:length(str)
    sum_densap(i)=sum(sum(ap_map(:,:,i)));
    sum_densba(i)=sum(sum(ba_map(:,:,i)));
end
% Setup A and Setup B
setups=[zeros(47,1);ones(100,1)];
% Slice orientation
slice_ori=[str(:).sliceOri];
%Fractions
frv=reshape([str(:).frac_vert],32,length(str))';
frh=reshape([str(:).frac_horz],32,length(str))';
L23fr=[sum(frv(:,3:5),2) sum(frv(:,19:21),2)];
L4fr=[sum(frv(:,6:7),2) sum(frv(:,22:23),2)];
L5fr=[sum(frv(:,8:11),2) sum(frv(:,24:27),2)];
%Differences of fractions
L5frt=L5fr;
L5frt(find(L5fr(:,1)==0),1)=NaN ;
L5frt(find(L5fr(:,2)==0),2)=NaN ;
diffL23fr=L23fr(:,1)-L23fr(:,2);
diffL4fr=L4fr(:,1)-L4fr(:,2);
diffL5fr=L5frt(:,1)-L5frt(:,2);
%span
[span(:,1),span(:,2),span(:,3),span(:,4),span(:,5),span(:,6)] = span_perLayer(ex_map,in_map,1:147)
%span=reshape([str(:).span],6,length(str))';
%Principal compponent scores
scores=reshape([str(:).PCs],6,length(str))';
%subpixel soma
somac=reshape([str(:).subpixel_soma],2,length(str))';
%Centroid measurements
% ang_exL23=reshape([str(:).ang_exL23],10,length(str))';
% ang_exL4=reshape([str(:).ang_exL4],10,length(str))';
% ang_exL5=reshape([str(:).ang_exL5],10,length(str))';
% ang_inL23=reshape([str(:).ang_inL23],10,length(str))';
% ang_inL4=reshape([str(:).ang_inL4],10,length(str))';
% ang_inL5=reshape([str(:).ang_inL5],10,length(str))';
%Pial depth
pia_input=[str(:).pialD]';
%Orientations
oris=[0:45:315];
%visual responses properties
%1:OSI,2:DSI,3:ODI,4:ORI,5:DIR,6:Capeak,7:Tuning Width, 8:SF,9
%TF,10:SAD,11:Noise,12:PCI
% temp=[str(:).ORIpref]+90
% temp(temp>179) = temp(temp>179) -180
% % for i=1:length(temp)
% %     if temp(i)>179.99==1
% %         oric(i)=temp(i)-180
% %     else
% %         oric(i)=temp(i);
% %     end
% % end
% oric=temp;
od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
             ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci]]';
%ipsi, contra, bino, unres
ipsi_id=find([str(:).ipsi]==1);
contra_id=find([str(:).contra]==1);
bino_id=find([str(:).bino]==1);
unres_id=find([str(:).unres]==1);
%% 
[ang_in] = centroid_map(in_map(:,:,:),somac(:,1),somac(:,2),[1:147],0);
[ang_inL23] = centroid_map(in_map(3:5,:,:),somac(:,1),somac(:,2),[1:147],2);
[ang_inL4] = centroid_map(in_map(6:7,:,:),somac(:,1),somac(:,2),[1:147],5);
[ang_inL5] = centroid_map(in_map(8:11,:,:),somac(:,1),somac(:,2),[1:147],7);
%EX ang centroid
[ang_exL23] = centroid_map(ex_map(3:5,:,:),somac(:,1),somac(:,2),[1:147],2);
[ang_exL4] = centroid_map(ex_map(6:7,:,:),somac(:,1),somac(:,2),[1:147],5);
[ang_exL5] = centroid_map(ex_map(8:11,:,:),somac(:,1),somac(:,2),[1:147],7);

%% Figure 1 panels
close all;
% Histogram of pia distribution with color coded in vivo morpho by itself
% for panel B
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 400, 400]);
histogram(pia_input,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
hold on;histogram(pia_input(morph_cells_id),'FaceColor','m','Orientation','horizontal');
hold on;histogram(pia_input(iviv_cells),'FaceColor','b','Orientation','horizontal');
legend([' Input (n=' num2str(length(pia_input)),')'],[' Morph. (n=' num2str(length(pia_input(morph_cells_id))),')']...
    ,[' In vivo (n=' num2str(length(pia_input(iviv_cells))),')']);
legend boxoff ;set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',12);
%Example input maps of cell 99 in Panel D
plot_maps(str,ex_map_raw,in_map_raw,[1:147],99,pia_input);
%Polar plot peak normalized in Panel C
u=99;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [500, 600, 200, 200]);
polarplot(deg2rad(oris([1:end 1])),str(u).TuningCurve([1:8 1])/max(str(u).TuningCurve));hold on;
polarplot(deg2rad(oris([1:end 1])),str(u).TuningCurve([9:end 1])/max(str(u).TuningCurve));
ax = gca;rticks([0:0.5:1]);ax.LineWidth = 2;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
 ax.ThetaTick = [0:45:360]; 
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure1\Final_panels\'
savepdf_SW(fn,1);
%% Figure 2 panels
close all;
% %Example input maps of cell 1,80,147 for panel A
%  plot_maps(str,ex_map_raw,in_map_raw,[1:147],1,pia_input);
%  plot_maps(str,ex_map_raw,in_map_raw,[1:147],80,pia_input);
%  plot_maps(str,ex_map_raw,in_map_raw,[1:147],147,pia_input);
% %Plot average maps all
% plot_avg_maps(str,1:147,ex_map,in_map,pia_input,10,0,[]);

%Fraction for EX and IN Vertical and Horizontal + differences for panel B
%and C: CALL FUNCTION display_inputs
[stats_g] = display_inputs([frv],[frh],frv(:,1:16)-frv(:,17:end),frh(:,1:16)-frh(:,17:end),[]);
%histograms for L23, L4, L5 Vertical and Horizontal (span), panel D and E
[statsout]=extent_diff(diffL23fr,diffL4fr,diffL5fr,span)
%Pial depth correlation with EX IN L23, L4, L5 for panel F
plot_fraction(L23fr,L4fr,L5fr,pia_input);

%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure2\Final_panels\'
savepdf_SW(fn,1);

%% Figure 3 panels
close all;
%Overview of centroid x any with respect to soma for EX and IN, Panel B
a=1:147;
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,0,[],[]);
%% 

%EX vs IN horizontal centroid CoMx for L2/3, L4, L5, Panel C
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 675]);
subplot(3,2,1);
plot(ang_exL23(:,3)*69-ang_exL23(:,1)*69,ang_inL23(:,3)*69-ang_inL23(:,1)*69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-100 100];
ref.YData=[-100 100];
ref.Color='k';box off;xlim([-100 100]);ylim([-100 100]);hold on;title('L23','FontWeight','normal');xticks([-100:50:100]);yticks([-100:50:100]);
subplot(3,2,3);
plot(ang_exL4(:,3)*69-ang_exL4(:,1)*69,ang_inL4(:,3)*69-ang_inL4(:,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];ylabel('IN C_{x} (µm)','Color','b')
ref.Color='k';box off;xlim([-150 150]);ylim([-150 150]);hold on;title('L4','FontWeight','normal');xticks([-150:75:200]);yticks([-150:75:150]);
subplot(3,2,5);
plot(ang_exL5(:,3)*69-ang_exL5(:,1)*69,ang_inL5(:,3)*69-ang_inL5(:,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;xlim([-300 300]);ylim([-300 300]);hold on;title('L5','FontWeight','normal');xticks([-300:150:300]);yticks([-300:150:300]);
xlabel('EX C_{x} (µm)','Color','r')
%EX vs IN verticalcentroid CoMY for L2/3, L4, L5, Panel C
subplot(3,2,2);
plot(ang_exL23(:,4)*-69-ang_exL23(:,2)*-69,ang_inL23(:,4)*-69-ang_inL23(:,2)*-69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];;xlim([-150 150]);ylim([-150 150]);xticks([-150:100:150]);yticks([-150:100:150]);
ref.Color='k';box off;;hold on;title('L23','FontWeight','normal');
subplot(3,2,4);
plot(((ang_exL4(:,4)*-69))-ang_exL4(:,2)*-69,((ang_inL4(:,4)*-69))-ang_inL4(:,2)*-69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');;xlim([-325 -25]);ylim([-325 -25]);
;ref= refline(1,0);
set(gca,'FontSize',10);ref.LineStyle='--';
ref.XData=[-325 -25];
ref.YData=[-325 -25];
ylabel('IN C_{y} (µm)','Color','b')
ref.Color='k';box off;hold on;title('L4','FontWeight','normal');xticks([-325:100:-25]);yticks([-325:100:-25]);
subplot(3,2,6);
plot(((ang_exL5(:,4)*-69))-ang_exL5(:,2)*-69,((ang_inL5(:,4)*-69))-ang_inL5(:,2)*-69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--';ref.XData=[-500 -50];
ref.YData=[-500 -50];xlim([-500 -50]);ylim([-500 -50]);xticks([-500:225:-50]);yticks([-500:225:-50]);

ref.Color='k';box off;hold on;title('L5','FontWeight','normal');
xlabel('EX C_{y} (µm)','Color','r');
%% Statisistcs
%Horizontal
[aa bb]=intersect(find(~isnan(ang_exL23(:,3)*69-ang_exL23(:,1)*69)),find(~isnan(ang_inL23(:,3)*69-ang_inL23(:,1)*69)));
[p1]=signrank(abs(ang_exL23(aa,3)*69-ang_exL23(aa,1)*69),abs(ang_inL23(aa,3)*69-ang_inL23(aa,1)*69))
[aa bb]=intersect(find(~isnan(ang_exL4(:,3)*69-ang_exL4(:,1)*69)),find(~isnan(ang_inL4(:,3)*69-ang_inL4(:,1)*69)));
[p2]=signrank(abs(ang_exL4(aa,3)*69-ang_exL4(aa,1)*69),abs(ang_inL4(aa,3)*69-ang_inL4(aa,1)*69))
[aa bb]=intersect(find(~isnan(ang_exL5(:,3)*69-ang_exL5(:,1)*69)),find(~isnan(ang_inL5(:,3)*69-ang_inL5(:,1)*69)));
[p3]=signrank(abs(ang_exL5(aa,3)*69-ang_exL5(aa,1)*69),abs(ang_inL5(aa,3)*69-ang_inL5(aa,1)*69))
%Vertical
[aa bb]=intersect(find(~isnan(ang_exL23(:,4)*69-ang_exL23(:,2)*69)),find(~isnan(ang_inL23(:,4)*69-ang_inL23(:,2)*69)));
[p4]=signrank(abs(ang_exL23(aa,4)*69-ang_exL23(aa,2)*69),abs(ang_inL23(aa,4)*69-ang_inL23(aa,2)*69))
[aa bb]=intersect(find(~isnan(ang_exL4(:,4)*69-ang_exL4(:,2)*69)),find(~isnan(ang_inL4(:,4)*69-ang_inL4(:,2)*69)));
[p5]=signrank(abs(ang_exL4(aa,4)*69-ang_exL4(aa,2)*69),abs(ang_inL4(aa,4)*69-ang_inL4(aa,2)*69))
[aa bb]=intersect(find(~isnan(ang_exL5(:,4)*69-ang_exL5(:,2)*69)),find(~isnan(ang_inL5(:,4)*69-ang_inL5(:,2)*69)));
[p6]=signrank(abs(ang_exL5(aa,4)*69-ang_exL5(aa,2)*69),abs(ang_inL5(aa,4)*69-ang_inL5(aa,2)*69))
%% 
figure;histogram(abs((ang_exL23(aa,4)*69-ang_exL23(aa,2)*69))-abs((ang_inL23(aa,4)*69-ang_inL23(aa,2)*69)))
%% 
[aa bb]=intersect(find(~isnan(ang_exL4(:,4)*69-ang_exL4(:,2)*69)),find(~isnan(ang_inL4(:,4)*69-ang_inL4(:,2)*69)));
%figure;histogram(abs((ang_exL4(aa,4)*69)-abs(ang_exL4(aa,2)*69)))-abs((ang_inL4(aa,4)*69)-abs(ang_inL4(aa,2)*69)))
%% 
aa=[];
[aa bb]=intersect(find(~isnan(ang_exL5(:,4)*69-ang_exL5(:,2)*69)),find(~isnan(ang_inL5(:,4)*69-ang_inL5(:,2)*69)));
figure;histogram(abs((ang_exL5(aa,4)*69-ang_exL5(aa,2)*69))-abs((ang_inL5(aa,4)*69-ang_inL5(aa,2)*69)))
%% 
aa=[]
[aa bb]=intersect(find(~isnan(ang_exL5(:,3)*69-ang_exL5(:,1)*69)),find(~isnan(ang_inL5(:,3)*69-ang_inL5(:,1)*69)));
figure;histogram(abs((ang_exL5(aa,3)*69-ang_exL5(aa,1)*69))-abs((ang_inL5(aa,3)*69-ang_inL5(aa,1)*69)))
%% 
figure;histogram(abs((ang_inL5(unres_id,4)*69-ang_inL5(unres_id,2)*69)));hold on;histogram(abs((ang_inL5(a,4)*69-ang_inL5(a,2)*69)))
%% 


%Principal component corr matrix
com=[];com=[scores]; 
O=correlation_matrix(com,0);title('');
fG=O(4:6,1:3)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 200])
imagesc(fG);c=colorbar;pos = get(c,'Position');[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:3]);yticks([1:1:3]);
xticklabels({'PC1in','PC2in','PC3in'});xtickangle(45);set(gca,'FontSize',10);xlabel('');
yticklabels({'PC1in','PC2in','PC3in'});ytickangle(45);set(gca,'FontSize',10);ylabel('');
%PC2in with Cy
corr_plot(scores(:,5),abs((ang_inL23(:,4)-ang_inL23(:,2))*-69),[],{'','',''});xlabel('PC2_{in}','Color','b');
%PC1ex with 
corr_plot(scores(:,1),pia_input,[],{'','',''});xlabel('PC1_{ex}','Color','r');
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);set(gca,'Ydir','reverse');ylim([100 400]);yticks([100:150:400]);xticks([-2:1:2]);
%CoMX vs PC13in, Panel G
corr_plot((ang_inL23(:,3)-ang_inL23(:,1))*69,scores(:,6),[],{'','',''});ylabel('PC3_{in}','Color','b');
xlabel('IN C_{x} L2/3 (µm)','Color','b');set(gca,'FontSize',10);ylim([-1.2 1.2]);yticks([-1.2:0.6:1.2]);xlim([-120 120]);xticks([-120:60:120])
%% Correlation matrix Pial depth Cx Cy layers for supplement
com=[];com=[ang_exL23(:,3)-ang_exL23(:,1)...
    ang_inL23(:,3)-ang_inL23(:,1) ang_exL23(:,4)-ang_exL23(:,2) ang_inL23(:,4)-ang_inL23(:,2)... 
    ang_exL4(:,3)-ang_exL4(:,1) ang_inL4(:,3)-ang_inL4(:,1) ang_exL4(:,4)-ang_exL23(:,2) ang_inL4(:,4)-ang_inL4(:,2)... 
    ang_exL5(:,3)-ang_exL5(:,1) ang_inL5(:,3)-ang_inL5(:,1) ang_exL5(:,4)-ang_exL5(:,2) ang_inL5(:,4)-ang_inL5(:,2)...
    pia_input]; 
O=correlation_matrix(com,0);title('');

xticks([1:1:13]);yticks([1:1:13]);caxis([-1 1]);c=colorbar;c.Ticks=[-1:0.5:1];
xticklabels({'C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','PiaD'});set(gca,'FontSize',10);xlabel('');set(gca,'FontSize',12)
yticklabels({'C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','PiaD'});set(gca,'FontSize',10);ylabel('');set(gca,'FontSize',12)
%% 
%Correlatoion matrix PCs and reald data 
com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1) L4fr(:,2)  ang_exL23(:,3)-ang_exL23(:,1)...
   ang_inL23(:,3)-ang_inL23(:,1) ang_exL23(:,4)-ang_exL23(:,2) ang_inL23(:,4)-ang_inL23(:,2) pia_input scores]
G=correlation_matrix(com,0);close(gcf);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);imagesc(G(10:15,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:18]);yticks([1:1:18]);
% xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});xtickangle(45)
% yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});ylabel('Feature');
yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}'});xtickangle(45);set(gca,'FontSize',12)
xticklabels({'L2/3fr','L4fr','L2/3fr','L4fr','C_{x} L2/3','C_{x} L2/3','C_{y} L2/3','C_{y} L2/3','Pial depth'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);
%% 
%Correlatoion matrix PCs and reald data 
com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1) L4fr(:,2)  abs(ang_exL23(:,3)-ang_exL23(:,1))...
   abs(ang_inL23(:,3)-ang_inL23(:,1)) abs(ang_exL23(:,4)-ang_exL23(:,2)) abs(ang_inL23(:,4)-ang_inL23(:,2)) pia_input scores]
G=correlation_matrix(com,0);close(gcf);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);imagesc(G(10:15,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:18]);yticks([1:1:18]);
% xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});xtickangle(45)
% yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});ylabel('Feature');
yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}'});xtickangle(45);set(gca,'FontSize',12)
xticklabels({'L2/3fr','L4fr','L2/3fr','L4fr','C_{x} L2/3','C_{x} L2/3','C_{y} L2/3','C_{y} L2/3','Pial depth'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12);
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure3\Final_panels\'
savepdf_SW(fn,1);
%% UMAP embedding
% get the pial depth
pialD=cat(1,str.pialD);
% get the layer 4 excitation
frac4ex = cat(1,str.frac_vert);
frac4ex = sum(frac4ex(:,6:7),2);
% get the inhibitory L23 angle

% get the inhibitory L23 x centroid
centroidX23in=abs(ang_inL23(:,3)-ang_inL23(:,1));
centroidY23in=abs(ang_inL23(:,4)-ang_inL23(:,2));
% assemble the feature vector
% cell_cell = cat(2,pialD,frac4ex,...
%     centroidX23in);
cell_cell = cat(2,pialD,frac4ex,...
     centroidX23in,centroidY23in);
% cell_cell = cat(2,pcs(:,2),ang);
cell_cell = normr_2(cell_cell,2);
%% Run UMAP on the data
% run the embedding from scratch
[reduced_data, umap] = run_umap(cell_cell, 'n_neighbors', 30, 'min_dist', 0.5);

%% Figure 4 panels
close all;
%UMAP projections Panel C and D
%load embedding
%load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\embedding_X_only\reduced_data.mat');
%in vitro values
plot_list = {'frac_vert','ang_inL23','pialD','span'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%in vivo values non cyclic 
plot_list = {'OSIpref','DSIpref','ODIpref','Sigmapref','Capeakpref'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%cyclic values
plot_list ={'ORIpref','DIRpref'}
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'phasemap');
autoArrangeFigures;

 %figure(1);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];
 figure(2);c=colorbar;caxis([0 0.7]);c.Ticks=[0:0.35:0.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(3);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(5);c=colorbar;caxis([0 110]);c.Ticks=[0:55:110];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(6);c=colorbar;caxis([0 160]);c.Ticks=[0:80:160];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(7);c=colorbar;caxis([120 360]);c.Ticks=[120:120:360];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(8);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% 
 figure(14);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(15);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(16);c=colorbar;caxis([-1 1]);c.Ticks=[-1:1:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(17);c=colorbar;caxis([2 4]);c.Ticks=[2:5];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(19);c=colorbar;caxis([0 180]);c.Ticks=[0:45:180];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(20);c=colorbar;caxis([0 360]);c.Ticks=[0:90:360];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% 
%  close(figure(1));close(figure(4));close(figure(9));close(figure(10));close(figure(11));close(figure(12));close(figure(13));
%  close(figure(8));close(figure(17));close(figure(18));
%% 
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200508_embedding\'
savepdf_SW(fn,1);
%% 
%Correlation Matrix Panel B
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) abs(ang_exL23(:,4)-ang_exL23(:,2))...
   abs(ang_inL23(:,4)-ang_inL23(:,2)) pia_input od_out_iviv(:,[1 2 3 4 5 6 7])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data OSI
a=find(od_out_iviv(:,1)>0.25)
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1) L4fr(a,2) abs(ang_exL23(a,3)-ang_exL23(a,1))  abs(ang_inL23(a,3)-ang_inL23(a,1))...
   abs(ang_exL23(a,4)-ang_exL23(a,2)) abs(ang_inL23(a,4)-ang_inL23(a,2)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
M=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data DSI
a=find(od_out_iviv(:,2)>0.25)
com=[];com=[L23fr(a,1) L23fr(a,2) L4fr(a,1)  L4fr(a,2) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) abs(ang_exL23(a,4)-ang_exL23(a,2)) abs(ang_inL23(a,4)-ang_inL23(a,2)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
O=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
close(gcf);
%Assemble matrix
fG=[G([10 11 12],1:9) ; M([13],1:9) ; O([14],1:9) ; G([15 16],1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 200])
imagesc(fG);c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]);
xticklabels({'L2/3_{ex}','L2/3_{in}','L4_{ex}','L4_{in}','C_{x} L2/3 ','C_{x} L2/3','C_{y} L2/3','C_{y} L2/3', 'Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI*','DIR*','TW','Ca_{peak}'});set(gca,'FontSize',10);
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
% c.Label.Position = [pos(1)/2 pos(2)+1]; % to change its position
% c.Label.Rotation = 0; % to rotate the text
%% 
%Selection of correlation from Panel B shown in Panel E
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),L4fr(a,1),[],{'','',''});ylabel('EX L4fr','Color','r');
xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0 180]);xticks([0:45:180]);
%% Cx and Cy with ori pref 
a=find(od_out_iviv(:,1)>0.25);
corr_plot((abs(ang_inL23(a,3)-ang_inL23(a,1)))*69,(abs(ang_inL23(a,4)-ang_inL23(a,2)))*69,od_out_iviv(a,4),{'','',''});ylabel('C_{y} (µm)','Color','b');xlabel('C_{x} (µm)','Color','b');
%xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0
%180]);xticks([0:45:180]);ylim([-10 150])
c=colorbar;c.Label.String = 'Orientation (deg)';caxis([0 180]);c.Ticks=[0:45:180];set(gca,'FontSize',10);
%% Small correlation matrix 
a=find(od_out_iviv(:,1)>0.25);
com=[];com=[abs(ang_inL23(a,3)-ang_inL23(a,1)) abs(ang_inL23(a,4)-ang_inL23(a,2)) pia_input(a) od_out_iviv(a,4)]; 
O=correlation_matrix(com,0);title('');xticks([1:1:4]);yticks([1:1:4]);
xticklabels({'C_{x}','C_{y}','Pial depth','ORI*'});xtickangle(45);set(gca,'FontSize',10);xlabel('');
yticklabels({'C_{x}','C_{y}','Pial depth','ORI*'});ytickangle(45);set(gca,'FontSize',10);ylabel('');
%% 
a=find(od_out_iviv(:,1)>0.25 & pia_input<300);
corr_plot(pia_input(w1),abs((ang_inL23(w1,3)-ang_inL23(w1,1))),od_out_iviv(w1,4),{'','',''});ylabel('C_{\alpha} (deg)','Color','b');

% %% Alternative polarplot
% a=find(od_out_iviv(:,1)>0.25);
% rho=abs(a_L23in(a));
% theta=od_out_iviv(a,4);%or theta=linspace(0,180,numel(rho));
% figure;set(gcf,'color','w')
% polarscatter(deg2rad(theta),rho,'ko','filled')
% ax = gca;ax.LineWidth = 1;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
% ax.ThetaLim = [0 180];
% ax.RLim = [0 150];
% % ax.RTickLabel = []; 
%  ax.ThetaTick = [0:45:180]; 
%  ax.RAxis.Label.String = 'My Label';
 %Examples for Fig 4A     %% 
%cell 126
id_m=126;
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp14783\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex1=Fit(2).contra.FittedDataOri
ind_tr1=peaks(2).contra;
fit_ex1oe=Fit(2).ipsi.FittedDataOri
ind_tr1oe=peaks(2).ipsi;
%cell 139
id_m2=139;
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp14875\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex2=Fit(3).ipsi.FittedDataOri;
ind_tr2=peaks(3).ipsi;
fit_ex2oe=Fit(3).contra.FittedDataOri
ind_tr2oe=peaks(3).contra;
shift_ori_oe(fit_ex1,fit_ex1oe,fit_ex2,fit_ex2oe,ind_tr1,ind_tr1oe,ind_tr2,ind_tr2oe,od_out_iviv(id_m,1),od_out_iviv(id_m2,1),od_out_iviv(id_m,7),od_out_iviv(id_m2,7));
plot_maps(str,ex_map_raw,in_map_raw,[1:147],139,pia_input);
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\Final_panels\'
savepdf_SW(fn,1);

%% Figure 5 panels

%Overview of centroid x and y with respect to soma for EX and IN and OSI colurcoded, Panel B
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
%% Vectro length vs ori pref
close all;
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),ang_exL23(a,8)*69,[],{'','',''});ylabel('CVL L2/3 (µm)');xlabel('Orientation (deg)');
yticks([0:50:100]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10);
corr_plot(od_out_iviv(a,4),ang_inL23(a,8)*69,[],{'','',''});ylabel('CVL L2/3 (µm)');xlabel('Orientation (deg)');
yticks([0:50:100]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10);
%% 
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),abs(ang_exL23(a,3)*69-ang_exL23(a,1)*69),[],{'','',''});

%% 
% ylabel('CVL L2/3 (µm)');xlabel('Orientation (deg)');
%yticks([0:50:100]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10);
corr_plot(od_out_iviv(a,4),ang_inL23(a,8)*69,[],{'','',''});ylabel('CVL L2/3 (µm)');xlabel('Orientation (deg)');
yticks([0:50:100]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10);

%% 
%Rolling average plots L23CVL
parameter_vector = abs(ang_exL23(:,8)*69)
parameter_vector2 = abs(ang_inL23(:,8)*69)
%parameter_vector = abs(ang_exL23(:,3)*69-ang_exL23(:,1)*69-(ang_inL23(:,3)*69-ang_inL23(:,1)*69))
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('CVL L2/3 (µm)','Color','k');
%   ylim([15 85]);yticks([15:35:85]);
%   xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
%% IN
a=find(od_out_iviv(:,1)>0.25);  
%par=abs(ang_inL4(a,3)*69-ang_inL4(a,1)*69)
par=ang_inL23(a,8)*69
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'15-65°','100-150°'});xtickangle(45);
ylabel('CVL L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);
%% EX
a=find(od_out_iviv(:,1)>0.25);  
%par=abs(ang_inL4(a,3)*69-ang_inL4(a,1)*69)
par=ang_exL23(a,8)*69
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'15-65°','100-150°'});xtickangle(45);
ylabel('CVL L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);

%% EX
a=find(od_out_iviv(:,1)>0.25);  
%par=abs(ang_inL4(a,3)*69-ang_inL4(a,1)*69)
par=ang_inL4(a,8)*69
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'15-65°','100-150°'});xtickangle(45);
ylabel('CVL L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);
%% 
data= vertcat(abs(ang_exL4(g2,4)*69-ang_exL4(g2,2)*69),abs(ang_inL4(g2,4)*69-ang_inL4(g2,2)*69));
groups_idx=vertcat(ones(length(g2),1)*1,ones(length(g2),1)*2)
[statsout] = barplot_sw(data,groups_idx,{'','Cvr L2/3'});ylim([0 100]);yticks([0:25:100])


%% Cx and Cy sperately

parameter_vector = abs((ang_inL23(:,3)-ang_inL23(:,1))*69)
parameter_vector2 = abs((ang_inL23(:,4)-ang_inL23(:,2))*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,65,2,1)
ylabel('IN Length (µm)','Color','b');
  %ylim([15 85]);yticks([15:35:85]);
  %xlim([0 180]);xticks([0:45:360]);
set(gca,'FontSize',10);
%% 
a=find(od_out_iviv(:,1)>0.25); 
par=abs(a,4)
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 
a=find(od_out_iviv(:,1)>0.25); 
par=max([span(a,[1 2 3]) span(a,[4 5 6])],[],2)
par=max([span(a,4)],[],2)
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 
parameter_vector = span(:,1)
parameter_vector2 = span(:,4)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('IN Length (µm)','Color','b');
  %ylim([15 85]);yticks([15:35:85]);
  %xlim([0 180]);xticks([0:45:360]);
set(gca,'FontSize',10);
%% 

par=ang_exL23(:,8)*69-ang_inL23(:,8)*69
%par=ang_inL23(:,8)*69
g1=find(od_out_iviv(:,5)>10 & od_out_iviv(:,5)<100) ;
g2=find(od_out_iviv(:,5)>180 & od_out_iviv(:,5)<270);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 
par=ang_exL4(:,8)
g1=find(od_out_iviv(:,5)>100 & od_out_iviv(:,5)<150) ;
g2=find(od_out_iviv(:,5)>270 & od_out_iviv(:,5)<320);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);

%% 
a=find(od_out_iviv(:,1)>0);  

corr_plot(ang_exL23(a,3)*69-ang_inL23(a,3)*69,od_out_iviv(a,5),od_out_iviv(a,4),{'','',''})
%corr_plot(ang_exL23(a,3)*69-ang_inL23(a,1)*69,od_out_iviv(a,2),od_out_iviv(a,3),{'','',''})
%% 
a=find(od_out_iviv(:,2)>0.25);  
figure;scatter(ang_exL23(a,3)*69,ang_inL23(a,3)*69,20,od_out_iviv(a,5),'filled');refline(1,0)
[cmap]=phasemap;colormap(cmap);axis equal;
%% 

parameter_vector = ang_exL23(:,8)*69-ang_inL23(:,8)*69

rolling_avg_display(str,parameter_vector,parameter_vector2,65,2,1)
ylabel('IN Length (µm)','Color','b');
%% GLM fitting


%% 
%a=find(od_out_iviv(:,1)>0.25 & abs(ang_inL4(:,3)*69-ang_inL4(:,1)*69)<300);
corr_plot(abs(ang_inL4(a,4)*69-ang_inL4(a,2)*69),od_out_iviv(a,4),pia_input(a),{'','',''})
a=find(od_out_iviv(:,1)>0.25);  
%par=abs(ang_inL4(a,3)*69-ang_inL4(a,1)*69)
par=ang_inL4(a,8)*69
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 



%% 

%Rolling average plots
parameter_vector = abs(ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,45,1);
ylabel('C_{x} L2/3 (µm)','Color','b');
ylim([5 60]);yticks([10:25:60]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
%L4 angle alpha 
parameter_vector = abs(a_L4in);
rolling_avg_display(str,parameter_vector,45,1);
ylabel('C_{\alpha} L4 (deg)','Color','b');
xlim([0 180]);xticks([0:45:180]);
ylim([2 20])
set(gca,'FontSize',10);

%dual barplots L23alpha
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=abs(a_L23ex(a));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45);
xlabel('EX','Color','r');
ylabel('C_{\alpha} L23 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 120]);yticks([0:60:120]);
%ylim([0 100]);
%IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(a_L23in(a));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');
ylabel('C_{\alpha} L23 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 120]);yticks([0:60:120]);
%dual barplots L23X
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_exL23(a,3)-ang_exL23(a,1))*69;
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('EX','Color','r');ylabel('C_{x} L23 (µm)','Color','k');set(gca,'FontSize',10); 
ylim([0 80]);

%IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL23(a,3)-ang_inL23(a,1))*69;
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');set(gca,'FontSize',10); 
xlabel('EX','Color','r');ylabel('C_{x} L23 (µm)','Color','k');set(gca,'FontSize',10); 
ylim([0 80]);;

%dual barplots L4alpha
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=abs(a_L4ex(a));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45);
xlabel('EX','Color','r');
ylabel('C_{\alpha} L4 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 40]);
%IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(a_L4in(a));;
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');
ylabel('C_{\alpha} L4 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 40]);

%dual barplots for L4 fraction in
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,1);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('L4fr','Color','r');set(gca,'FontSize',10);xtickangle(45)
%dual barplots for L4 fraction EX
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,2);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
;ylabel('L4fr','Color','b');set(gca,'FontSize',10);xtickangle(45)
%pia depth
a=find(od_out_iviv(:,1)>0.25);  
par=pia_input(a);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45)
%Correlation L4fr and CoMx/CoMalpha
% a=find(od_out_iviv(:,1)>0.25); 
% corr_plot(90-abs(ang_inL23(a,5)),L4fr(a,1),[],{'','',''});xlabel('L23 CoM\alpha (deg)','Color','b');
% ylabel('L4fr','Color','r');set(gca,'FontSize',10);
%% 
 a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL4(a,3)-ang_inL4(a,1))
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45)
%% Plot average maps all
w1=[];w2=[];
w1=find(od_out_iviv(:,4)>15 & od_out_iviv(:,4)<65 & od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);
w1a=find(od_out_iviv(:,5)>195 & od_out_iviv(:,5)<245 & od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25);
w1b=find(od_out_iviv(:,5)>15 & od_out_iviv(:,5)<65 & od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25);
w2=find(od_out_iviv(:,4)>100 & od_out_iviv(:,4)<150 & od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);
%  w1=find(pia_input<220);
%  w2=find(pia_input>220);
%% responsive or not
nvr=find([str(:).sftf_resp]==0 & [str(:).resp]==0)
nvra=find([str(:).resp]==0);
nvrb=find([str(:).sftf_resp]==1 & [str(:).resp]==0);
%nvr=find([str(:).resp]==0);
vr=find([str(:).resp]==1);

%% 

a=find(od_out_iviv(:,1)>0.25);
corr_plot((abs(ang_inL23(a,3)-ang_inL23(a,1)))*69,(abs(ang_inL23(a,4)-ang_inL23(a,2)))*69,od_out_iviv(a,4),{'','',''});ylabel('C_{y} (µm)','Color','b');xlabel('C_{x} (µm)','Color','b');
%xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0
%180]);xticks([0:45:180]);ylim([-10 150])
c=colorbar;c.Label.String = 'Orientation (deg)';caxis([0 180]);c.Ticks=[0:45:180];set(gca,'FontSize',10);
%% 

figure;histogram((ang_inL23(nvr,8))*69,8);hold on;histogram((ang_inL23(w2,8))*69,8)
%% 
figure;histogram(od_out_iviv(w1a,5))

%% 

corr_plot(ang_inL23(nvr,3)-ang_inL23(nvr,1),ang_inL23(nvr,4)-ang_inL23(nvr,2),[],{'','',''});ylabel('C_{y} (µm)','Color','b');xlabel('C_{x} (µm)','Color','b');
%xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0
%180]);xticks([0:45:180]);ylim([-10 150])


%% 

 %par=ang_exL23(:,4)-ang_exL23(:,2)
 %par=ang_inL23(:,8)
 %par=ang_inL23(:,8)
%par=reshape(sum(sum(in_map(4:5,:,:))),1,147)
%par=span(:,4)*69
par=ang_exL23(:,8)
g1=nvr
g2=vr
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'NV','VR'});
ylabel('Span','Color','k');set(gca,'FontSize',10);xtickangle(45)
%% 
par=od_out_iviv(:,13)
g1=w1
g2=w2
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'non visual resp','100-150°'});
ylabel('Cvr (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45)
%% 

a=find(od_out_iviv(:,7)>28);  
par=ang_exL23(:,8)*69
g1=a
g2=w2
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'broad tuning width','100-150°'});
ylabel('Vector length  (µm)','Color','r');set(gca,'FontSize',10);xtickangle(45)
%% 
s_vr=[];
inpuj=ang_exL23
tz1=find(inpuj(:,10)==1 | inpuj(:,10)==3);
tz2=find(inpuj(:,10)==2 | inpuj(:,10)==4);
 para1=inpuj(tz1,8)*-1;
para2=inpuj(tz2,8)*1;
s_vr=ones(147,1)*NaN;
s_vr(tz1)=para1;
s_vr(tz2)=para2;

%% Direction
%Rolling average plots L23alpha
%parameter_vector = s_vr*69
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69
rolling_avg_display(str,parameter_vector,parameter_vector2,65,2,0);title('');
 ylabel('C_{x} L2/3 (µm)','Color','r');xlabel('Direction (deg)');
 %ylim([-70 100]);%yticks([-70:70:140]);
 xlim([0 360]);xticks([0:90:360]);
 set(gca,'FontSize',10);
 %% 
a=find(od_out_iviv(:,2)>0.25 & od_out_iviv(:,2)>0);  
%par=s_vr*69
par=ang_exL23(a,3)-ang_exL23(a,1)
g1=find(od_out_iviv(a,5)>75 & od_out_iviv(a,5)<140)
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<320)
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~125','~305°'});
ylabel('C_{x} L2/3 (µm)','Color','r');set(gca,'FontSize',10);xtickangle(45)
%% 

data= vertcat(ang_inL23(nvr,8)*69,ang_inL23(w1,8)*69,ang_inL23(w2,8)*69);
groups_idx=vertcat(ones(length(nvr),1),ones(length(w1),1)*2,ones(length(w2),1)*3)
[statsout] = barplot_sw(data,groups_idx,{'','Cvr L2/3'});ylim([0 100]);yticks([0:25:100])

%% 
par=ang_inL23(:,8)*69
g1=contra_id
g2=ipsi_id
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 
wz=find(od_out_iviv(:,4)>75 & od_out_iviv(:,1)>0.25)
corr_plot(od_out_iviv(wz,7),(ang_inL23(wz,8))*69,[],{'Tuning Width','Cvr',''}); 
%% 
figure;scatter(abs((ang_inL23(w1,3)-ang_inL23(w1,1))*69),abs((ang_inL23(w1,4)-ang_inL23(w1,2))*69),20,od_out_iviv(w1,7),'filled');hold on;xlim([0 80]);ylim([0 80])
figure;scatter(abs((ang_inL23(w2,3)-ang_inL23(w2,1))*69),abs((ang_inL23(w2,4)-ang_inL23(w2,2))*69),20,od_out_iviv(w2,7),'filled');hold on;xlim([0 80]);ylim([0 80])
%% 

corr_plot(ang_inL23([w1; w2],8),od_out_iviv([w1 ;w2],4),pia_input([w1 ;w2]),{'','',''}); 
%% 
for j=1:length(in_map)
    if ang_inL23(j,3)-ang_inL23(j,1)<0
f_inmap(:,:,j)=fliplr(in_map(:,:,j));
frh_f(j,:)=fliplr(frh(j,:));
    else
        f_inmap(:,:,j)=in_map(:,:,j);
        frh_f(j,:)=frh(j,:);
    end
end
%% 
figure;imagesc(nanmean(in_map(3:5,:,w2),3))
%% 
figure;imagesc(nanmean(f_inmap(3:5,:,w1),3))
%% 
t1=nanmean(sum(f_inmap(3:5,:,w1)),3);
t2=nanmean(sum(f_inmap(3:5,:,w2)),3);


%% 
figure;plot(t1,'-r');hold on;plot(t2,'-k')
legend('~45','~135')
set(gcf,'color','w')
ylabel('Sum fraction L23');
xlabel('Pixel')
%% 

figure;plot(nanmean(frh_f(w2,16:end)',2),'-k');hold on;plot(nanmean(frh_f(w1,16:end)',2),'-r')
legend('~45','~135')
set(gcf,'color','w')
ylabel('Sum fraction all map');
xlabel('Pixel')
%% 
figure;plot(nanmean(frv(w1,16:end),1)');hold on;plot(nanmean(frv(w2,16:end),1)');
%% 


idx_w=NaN*ones(147,1);
idx_w(w1)=1;idx_w(w2)=2;
plot_avg_maps(str,1:147,ex_map,in_map,pia_input,1,1,idx_w,ang_inL23);
%% 
a=find(od_out_iviv(:,1)>0.25); 
corr_plot(od_out_iviv(a,5),(ang_inL23(a,3)-ang_inL23(a,1))*69,[],{'','',''}); 
%% 
corr_plot(od_out_iviv(w2,2),(abs(L23_t(3,w2)'-ang_inL23(w2,1))*69),[],{'','',''}); 
%% 
corr_plot(abs((ang_inL23(a,3)-ang_inL23(a,1))),od_out_iviv(a,4),[],{'','',''}); 
%% 
corr_plot(od_out_iviv(a,4),abs(ang_inL23(a,5)),[],{'','',''}); 
%% 
corr_plot(abs((ang_inL23(a,3)-ang_inL23(a,1))),pia_input(a),[],{'','',''}); 
%% 
figure;scatter3(abs((ang_inL23(a,3)-ang_inL23(a,1))),abs((ang_inL23(a,4)-ang_inL23(a,2))),od_out_iviv(a,4))
%% 
corr_plot(abs(ang_inL23(a,5)),od_out_iviv(a,4),[],{'','',''});
%% 

figure;
scatter(ang_inL23(w1,3)-ang_inL23(w1,1),ang_inL23(w1,4)-ang_inL23(w1,2),'k*')
hold on;
scatter(ang_inL23(w2,3)-ang_inL23(w2,1),ang_inL23(w2,4)-ang_inL23(w2,2),'g*');
xlim([-1 1])
%% 
a=find(od_out_iviv(:,1)>0.25);  
parameter_vector = abs(scores(:,6))
%parameter_vector = pia_input
rolling_avg_display(str,parameter_vector,45,1);title('');
 ylabel('|C_{\alpha} L2/3 (deg)|','Color','b');xlabel('Direction (deg)');
 
%% Direction
%Rolling average plots L23alpha
parameter_vector = a_L23in;
rolling_avg_display(str,parameter_vector,65,2);title('');
 ylabel('|C_{\alpha} L2/3 (deg)|','Color','b');xlabel('Direction (deg)');
 ylim([-70 120]);yticks([-70:70:140]);
 xlim([0 360]);xticks([0:90:360]);
 set(gca,'FontSize',10);
parameter_vector = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,65,2);
ylabel('|C_{x} L2/3 (µm)|','Color','b');xlabel('Direction (deg)');
title('');
 ylim([-60 70]);yticks([-60:60:70]);
 xlim([0 360]);xticks([0:90:360]);
 set(gca,'FontSize',10);
%% Direction barplot
a=find(od_out_iviv(:,2)>0.25);  
par=(ang_inL23(a,4)-ang_inL23(a,2))*69;
%par=a_L23in(a)
g1=find(od_out_iviv(a,5)>80 & od_out_iviv(a,5)<145) ;
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<325);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~135°','~315°'});
ylabel('C_{x} L2/3 (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45);
ylim([-100 100]);yticks([-100:50:100]);
%% 
a=find(od_out_iviv(:,2)>0.25);  
par=pia_input(a);
%par=a_L23in(a)
g1=find(od_out_iviv(a,5)>80 & od_out_iviv(a,5)<145) ;
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<325);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~135°','~315°'});
ylabel('C_{x} L2/3 (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45);

%% 

a=find(od_out_iviv(:,2)>0.25);  
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%par=a_L23in(a)
g1=find(od_out_iviv(a,5)>80 & od_out_iviv(a,5)<145) ;
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<325);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~135°','~315°'});
ylabel('C_{x} L2/3 (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45);
ylim([-100 100]);yticks([-100:50:100]);
%% 
a=find(od_out_iviv(:,2)>0.25);  
par=a_L23in(a)
g1=find(od_out_iviv(a,5)>80 & od_out_iviv(a,5)<145) ;
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<325);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~135°','~315°'});
ylabel('C_{alpha} L2/3 (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45);
ylim([-100 120]);yticks([-100:50:100]);

a=find(od_out_iviv(:,2)>0.25);  
par=a_L23ex(a)
g1=find(od_out_iviv(a,5)>80 & od_out_iviv(a,5)<145) ;
g2=find(od_out_iviv(a,5)>255 & od_out_iviv(a,5)<325);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'~135°','~315°'});
ylabel('C_{alpha} L2/3 (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45);
ylim([-100 120]);yticks([-100:50:100]);
%% 
a=find(od_out_iviv(:,2)>0.25); 
corr_plot(od_out_iviv(a,5),(ang_inL23(a,3)-ang_inL23(a,1))*69,[],{'','',''}); 
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure5\Final_panels\'
savepdf_SW(fn,1);
%% Perform sholl anaylsis, this should be incoorporated into the structure
%Get all morphtraces
for i=1:length(str)
    if ~isempty(str(i).morph)==1 
    zz{:,i}=str(i).morphtraces;
    else
        zz{:,i}=NaN;
    end
end
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147,0);
%% 
for i=1:length(zz)
if morph_cells(i)==1
oo=stats_tree(zz{1, i}{1, 1},{'HSE','HSN','VS2','VS3','VS4'},[],'-w -x');
op=stats_tree(zz{1, i}{1, 2},{'HSE','HSN','VS2','VS3','VS4'},[],'-w -x');
bnr(i)=length(oo.dstats.blen{:})
bnr_ba(i)=length(op.dstats.blen{:})
else
 bnr(i)=NaN;   
 bnr_ba(i)=NaN;   
end

end
%% Figure 6 panels
close all;
%Plot representative examples, cell 143 and 122 for panel A, 
%NOTE tuning curve: individual traces not in the str only the average atm
id_m=143;
plot_morphologies(str,id_m,1,1,1);
id_m2=122;
plot_morphologies(str,id_m2,2,1,1);
%Plot the tuning curves and width from the folders
%cell 143
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp15006\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex1=Fit(1).ipsi.FittedDataOri
ind_tr1=peaks(1).ipsi;
%cell 122
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\exp14757\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex2=Fit(2).contra.FittedDataOri;
ind_tr2=peaks(2).contra;
%call function shift_ori for plotting PANEL A,B
shift_ori(fit_ex1,fit_ex2,ind_tr1,ind_tr2,od_out_iviv(id_m,1),od_out_iviv(id_m2,1),od_out_iviv(id_m,7),od_out_iviv(id_m2,7));
%Plot the example sholl analyis for cell 143 and 122, PANEL C
temp=zz{143};
figure;
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(temp{1,1}, 20, '-s');
[sb, ddb, sdb, XP, YP, ZP, iD] = sholl_tree(temp{1,2}, 20, '-s');
temp=zz{122}
[s1, dd1, sd1, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,1}, 20, '-s');
[s1b, dd1b, sd1b, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,2}, 20, '-s');
close(gcf);

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 500, 250, 225]);
;plot(dd,s,'-k');box off;
hold on;plot(ddb,sb,'--k')
hold on;plot(dd1,s1,'-b');
hold on;plot(dd1b,s1b,'--b');
legend('Apical','Basal');legend boxoff;
ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%Plot correlations
%subsample cells with morph and fct. in vivo
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1 & str(i).resp==1 & ~isnan(str(i).OSIpref)==1
        m_res(i)=1;     
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=find(m_res==1);
%plot apical nr of brainch points vs OSI, PANEL D
corr_plot(morph_parameters(morph_res_sub,4),od_out_iviv(morph_res_sub,1),[],{'','',''});xlim([0 25]);ylim([0 1]);yticks([0:0.2:1]);
ylabel('OSI','Color','k');xlabel('Nr of branch points','Color','k');set(gca,'FontSize',10);
%plot peak nr of sholl crossings APICAL vs TW, PANEL D
%remove obvious TW outlier
morph_res_sub2=morph_res_sub;
morph_res_sub2(find(od_out_iviv(morph_res_sub2,7)>50))=[];
%Plot
corr_plot(max_s(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 15]);
ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
corr_plot(morph_parameters(morph_res_sub2,14),od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 15]);
ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',10)

%plot peak nr of sholl crossings BASAL vs TW, PANEL E
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 30]);
ylabel('Tuning Width','Color','k');xlabel('Peak Nr. sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%% 

%Plot correlation matrix, PANEL F
df=[morph_parameters(morph_res_sub,9) morph_parameters(morph_res_sub,10)];
db=[morph_parameters(morph_res_sub,19) morph_parameters(morph_res_sub,20)];
com=[];com=[morph_parameters(morph_res_sub,2) nanmax(df,[],2) dis_s(morph_res_sub)'  morph_parameters(morph_res_sub,4)  max_s(morph_res_sub)' ...
    morph_parameters(morph_res_sub,12) nanmax(db,[],2) dis_s_ba(morph_res_sub)'  morph_parameters(morph_res_sub,14) max_s_ba(morph_res_sub)'  od_out_iviv(morph_res_sub,[1 2 3 4 5 7 8]) pia_input(morph_res_sub)]
G=correlation_matrix(com,0);close(gcf)
Gf=G(11:end,1:10)
Gf(6,1)=0;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(Gf);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:8]);
 
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%% 
% Plot the input vs morphology density,panel G
for i=1:length(str)
%whole map
tmp1=ex_map(3:end,:,i);
tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:));sum(tmp2(:))];
end
corr_plot(max_s(morph_cells_id)',tot_input(morph_cells_id,2),[],{'','',''});xlabel('Peak Nr. sholl crossings','Color','k');
ylabel('Total input','Color','r');
%% Sum input for ex and in
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [900, 500, 200, 200])
scatter(max_s(morph_cells_id)',tot_input(morph_cells_id,2),'b.');
hold on;scatter(max_s(morph_cells_id)',tot_input(morph_cells_id,1),'r.');
xlabel('Peak nr. Sholl crossings','Color','k');
ylabel('Total input');set(gca,'FontSize',10);
%% Histograms for ranges
figure;histogram(od_out_iviv(morph_res_sub,4));
figure;histogram(od_out_iviv(morph_res_sub,1));
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure6\Final_panels\'
savepdf_SW(fn,1);
%% 



%% 







%% Supplementary Figures
%in vivo in vitro comparison
% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 [od_out sftf_out sftf_sel sftf_pref spon_out pia_all] = concat_invivo(L23_PC);
%% Supplementary Fig 3 Panel A
close all
%OSI and ORI
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 600, 200]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(od_out(find(~isnan(od_out(:,1))),1),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,1))),1),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('OSI','FontSize',12);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
xlim([-0.1 1]);xticks([0:0.2:0.8])
a=od_out_iviv(:,1)>0.25;
b=od_out(:,1)>0.25;set(gca,'FontSize',10);
subplot(1,2,2);
binRange = 0:45:179;
hcx = histcounts(od_out(b,4),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,4),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Orientation (deg)','FontSize',12);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);set(gca,'FontSize',10);
%DSI and DIR
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [800, 200, 600, 200]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(od_out(find(~isnan(od_out(:,2))),2),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,2))),2),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('DSI','FontSize',12);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);xlim([-0.1 1]);xticks([0:0.2:0.8])
a=od_out_iviv(:,2)>0.25;
b=od_out(:,2)>0.25;set(gca,'FontSize',10);
subplot(1,2,2);
binRange = 0:45:315;
hcx = histcounts(od_out(b,5),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,5),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Direction (deg)','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
xtickangle(45)
ylim([0 0.4]);yticks([0:0.2:0.4]);set(gca,'FontSize',10); 
%ipsi contra bino
for i=1:length(od_out)
if od_out(i,9)==1 & od_out(i,10)==0;
    o_all(i)=1;
elseif od_out(i,9)==0 & od_out(i,10)==1;
    o_all(i)=2;
elseif od_out(i,9)==1 & od_out(i,10)==1;
    o_all(i)=3;
else od_out(i,9)==0 & od_out(i,10)==0;
    o_all(i)=0;
end
end
fig22=figure;set(gcf,'color','w');set(fig22, 'Position', [300, 700, 250, 200]);set(gcf,'color','w');
resp_all=length(find(o_all==1))+length(find(o_all==2))+length(find(o_all==3));
resp_iviv=length(contra_id)+length(ipsi_id)+length(bino_id);
b=bar([1 2],[length(find(o_all==1))/resp_all length(find(o_all==2))/resp_all length(find(o_all==3))/resp_all; length(contra_id)/resp_iviv ...
    length(ipsi_id)/resp_iviv length(bino_id)/resp_iviv],'stacked');box off;ylabel('Fraction of cells');
 b(1, 1).FaceColor='b';b(1, 2).FaceColor='r';b(1, 3).FaceColor='w';
 xticklabels({'all','iviv'});legend('contra','ipsi','bino');legend boxoff;set(gca,'FontSize',10);
 %ODI distribution
% fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 600, 200]);set(gcf,'color','w');
% binRange =-1:0.2857:1;
% hcx = histcounts(od_out(find(~isnan(od_out(:,3))),3),[binRange Inf],'Normalization','probability');
% hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,3))),3),[binRange Inf],'Normalization','probability');
% b1=bar(binRange,[hcx;hcy]');box off;xlabel('ODI','FontSize',10);ylabel('Fraction of cells');
% b1(1).FaceColor=[0 0 0];
% b1(2).FaceColor=[0 0.7 0.6];
% ylim([0 0.4]);yticks([0:0.2:0.6]);set(gca,'FontSize',10);
 %ODI distribution
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 300, 200]);set(gcf,'color','w');
h=histogram(od_out(find(~isnan(od_out(:,3))),3),7,'Normalization','probability');h.FaceColor=[0 0 0];
hold on;h=histogram(od_out_iviv(find(~isnan(od_out_iviv(:,3))),3),7,'Normalization','probability');
h.FaceColor=[0 0.7 0.6];box off;
ylim([0 0.3]);yticks([0:0.15:0.3]);xlabel('ODI');ylabel('Fraction of cells');set(gca,'FontSize',10);

%TW distribution 
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 250, 200]);set(gcf,'color','w');
binRange = 0:20:90;
hcx = histcounts(od_out(find(~isnan(od_out(:,7))),7),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,7))),7),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Tuning width (deg)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.6]);yticks([0:0.2:0.6]);set(gca,'FontSize',10);

%Ca peak
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 250, 200]);set(gcf,'color','w');
binRange = 0:50:400;
hcx = histcounts(od_out(find(~isnan(od_out(:,6))),6),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,6))),6),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Ca_{peak} (?R/R0)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.7]);yticks([0:0.35:0.7]);set(gca,'FontSize',10);xtickangle(45);


%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp3\Final_panels\'
savepdf_SW(fn,1);
%% SPon activities display
%Spon activity
fig23 = figure;
set(fig23, 'Name', 'Binocular protocol');set(fig23, 'Position', [1000, 800, 600, 200]);set(gcf,'color','w');
subplot(1,3,1);
binRange = 0.0013:0.025:0.2;
hcx = histcounts(spon_out(find(~isnan(spon_out(:,1))),1),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,10))),10),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('SAD (Hz)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
xticks([0:0.05:0.15]);xtickangle(45);set(gca,'FontSize',10);
ylim([0 0.8]);yticks([0:0.4:0.8]);

subplot(1,3,2);
binRange = -7:1:2;
hcy = histcounts(log(abs(od_out_iviv(:,12))),[binRange Inf],'Normalization','probability');
hcx = histcounts(log(abs(spon_out(:,2))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('PCI','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
xticks([-7:2:2]);xtickangle(45);
ylim([0 0.8]);yticks([0:0.4:0.8]);set(gca,'FontSize',10);

%noise
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\noise\noise_correlations.mat')
noise_OD(2115:2115+39,:,:)=[];
for i=1:length(noise_OD)
    if o_all(i)==0
        noise_all(i)=NaN;
    elseif o_all(i)==1
        noise_all(i)=nanmean(noise_OD(i,:,1));
         elseif o_all(i)==2
        noise_all(i)=nanmean(noise_OD(i,:,2));
    else o_all(i)==3
        if od_out(i,3)>0
        noise_all(i)=nanmean(noise_OD(i,:,1));
        else od_out(i,3)<0
            noise_all(i)=nanmean(noise_OD(i,:,2));
        end
    end    
end

subplot(1,3,3);
binRange = -0.02:0.02:0.16;
hcx = histcounts(noise_all(find(~isnan(noise_all))),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,11))),11),[binRange Inf],'Normalization','probability');

b1=bar(binRange,[hcx;hcy]');box off;xlabel('Noise','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];xticks([-0.02:0.04:0.16]);xtickangle(45);
ylim([0 0.8]);yticks([0:0.4:0.8]);set(gca,'FontSize',10);
%% Corrrltaion with dark and input maps

%Correlation matirx PCI, SAD, noise with input
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) abs(ang_exL23(:,5)) abs(ang_inL23(:,5)) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) pia_input od_out_iviv(:,[10 11 12])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);

%Assemble matrix
fG=[G([10 11 12],1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 200])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]);
xticklabels({'L2/3_{ex}','L2/3_{in}','L4_{ex}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'SAD','noise','PCI'});set(gca,'FontSize',10);
c.Label.String = 'R'; c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10)

%corr_plot(abs(ang_exL23(:,3)-ang_exL23(:,1)),od_out_iviv(:,10),[],{'','',''});
%% Morphology and input maps
close all;
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) abs(ang_exL23(:,5)) abs(ang_inL23(:,5)) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) pia_input morph_parameters(:,1:21)]% morph_parameters(:,1) morph_parameters(:,8) dis_s'  morph_parameters(:,4) max_s'   nanmax(df,[],2)...
    %morph_parameters(:,12) morph_parameters(:,11) morph_parameters(:,18) dis_s_ba'  morph_parameters(:,14) max_s_ba'  nanmax(db,[],2)]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%Assemble matrix
fG=[G(10:end,1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 800, 600])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:16]);yticks([1:1:25]);
xticklabels({'L2/3_{ex}','L2/3_{in}','L4_{ex}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'RDAmax','LAtotal','PLAmax','BPA','BOAmax','BLA',...
                'PLA','WHA','XSA','YSA','RDBmax','LBtotal','PLBmax','BPB',...
                'BOBmax','BLB','PLB','WHB','XSB','YSB','NB'});set(gca,'FontSize',10);
c.Label.String = 'R'; c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10)
%% 
corr_plot(span(:,5)',morph_parameters(:,19) ,[],{'','',''});
%% 
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) 90-abs(ang_exL23(:,5)) 90-abs(ang_inL23(:,5)) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) pia_input morph_parameters(:,2) nanmax(df,[],2) dis_s'  morph_parameters(:,4)  max_s'  sum_densap' ...
    morph_parameters(:,12) nanmax(db,[],2) dis_s_ba'  morph_parameters(:,14) max_s_ba'  sum_densba']
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
fG=[G(10:end,1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 800, 600])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:16]);yticks([1:1:25]);
xticklabels({'L2/3_{ex}','L2/3_{in}','L4_{ex}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}'});xtickangle(45);set(gca,'FontSize',12)

corr_plot(L4fr(:,1),nanmax(db,[],2),[],{'','',''});
%% 
temp=zz{3}
[s1, dd1, sd1, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,1}, 20, '-s');
[s1b, dd1b, sd1b, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,2}, 20, '-s');
close(gcf);

%% Plot sholl crossings
close all;
figure(30);set(gcf,'color','w');set(figure(30), 'Position', [200, 200, 300, 200])
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147,1)
figure(30);
hold on;plot(dd1,s1,'-k');
hold on;plot(dd1b,s1b,'-m');
legend('Apical','Basal');legend boxoff;
ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
ylabel('Nr. of crossings');xlabel('Distance from Soma (µm)');
%% Example cell with sholl rings
id_m=3;
plot_morphologies(str,id_m,1,1,1);
ml=13
  xlim([-400 400]);
 ylim([-200 600]);
for i=1:ml
hold on;
  viscircles([0 pia_input(id_m)],20*i,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
 axis off
%% Histograms of parameters for morpho
morph_sel=[morph_parameters(:,[1 2]) max_s' dis_s' morph_parameters(:,[8]) morph_parameters(:,[11 12]) max_s_ba' dis_s_ba' morph_parameters(:,[18])]

 stri={'RDA_{max}(µm)','Total Length (µm)','Peak Nr. crossing','Dis. peak branch (µm)','WHA',...
     'RDB_{max}(µm)','Total Length (µm)','Peak Nr. crossing','Dis. peak branch(µm)','WHA'}
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200])
 for i=1:10
hold on;
subplot(2,5,i)
h=histogram(morph_sel(:,i),8,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5)
h.EdgeColor = 'k';
h.FaceColor = [0.5 0.5 0.5];
xlabel(stri(i));
ylim([0 50]);
xlim([0 max(morph_sel(:,i))+max(morph_sel(:,i))*0.25]);
hAxis = gca;
hAxis.YAxisLocation = 'left';    % 'left' (default) or 'right'
hAxis.XAxisLocation = 'bottom'
box off
%ylabel('Counts');
 end
 %% Example morphologies
 %mid=3, L3:7; L2:18
id_m=18;
plot_morphologies(str,id_m,1,1,1);
id_m=3;
plot_morphologies(str,id_m,1,1,1);
id_m=7;
plot_morphologies(str,id_m,1,1,1);
%% Correlation matrix between input, morpho, pia

df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2)  span(:,1) span(:,4) tot_input(:,1) tot_input(:,2) pia_input morph_sel]
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
fG=[G(10:end,1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 400])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);c.Ticks=[-1:0.5:1];
xticks([1:1:16]);yticks([1:1:25]);
xticklabels({'L2/3fr','L2/3fr','L4fr','L4fr','Span L2/3','Span L2/3', 'Total norm. input','Total norm. input','Pial depth'});xtickangle(45);set(gca,'FontSize',10)
yticklabels({'Radial dis.','Total Length','Peak Nr. crossing','Dis. peak branch','Width/Height',...
     'Radial dis.','Total Length','Peak Nr. crossing','Dis. peak branch','Width/Height'});set(gca,'FontSize',10)
%% Plot 2 significant correlations
a=find(~isnan(morph_sel(:,5)))
corr_plot(morph_sel(a,9),max(span(a,4:6)*69,[],2),[],{'','',''});ylabel('Max Span (µm)');xlabel('Dis. peak branch(µm)')
%corr_plot(nanmax(db(a,:),[],2),max(span(a,4:6),[],2),[],{'','',''});
corr_plot(morph_sel(a,5),pia_input(a),L4fr(a,1),{'','','EX L4fr'});ylabel('Pial depth (µm)');
xlabel('Width/Height');ylim([130 350]);yticks([130:100:350])
c=colorbar;caxis([0 0.8]);c.Ticks=[0:0.2:0.8];set(gca,'Ydir','reverse');
%% 
a=find(~isnan(morph_sel(:,5)))
corr_plot(morph_sel(a,10),max(span(a,1)*69,[],2),[],{'','','EX L4fr'});

%% Cell plotter, morphology, in vivo, maps
close all;
for i=67:69;
figure;
iviv_plotter(str,i)
end
%% Save individual morpho, in vivo, maps panels
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp1\Block1\'
savepdf_SW(fn,0);
%% Example for with and without TTX
pathName='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\TTX\161116\SW0005_inhibition_before\map04'
load_raw_map(pathName,str,68)
%% PCA supplementary
close all;
map_align_PCA(str);
%% PC score display sorted
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);
subplot(1,2,1)
[so soi]=sort(scores(:,1),'descend');
imagesc(scores(soi,1:3));c=colorbar;caxis([-2 2]);c.Ticks=[-2:1:2];
xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}'});xtickangle(45);ylabel('Cells');ylim([1 147]);
yticks([1:24:147]);set(gca,'FontSize',10);
subplot(1,2,2)
[so soi]=sort(scores(:,4),'descend');
imagesc(scores(soi,4:6));c=colorbar;caxis([-2 2]);c.Ticks=[-2:1:2];
xticklabels({'PC1_{in}','PC2_{in}','PC3_{in}'});xtickangle(45);ylim([1 147]);
yticks([1:24:147]);set(gca,'FontSize',10);
%% Display correlation between all PCs Multiple comparison 
com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1) L4fr(:,2) abs(a_L23ex) abs(a_L23in) ang_exL23(:,3)-ang_exL23(:,1)...
   ang_inL23(:,3)-ang_inL23(:,1) pia_input scores]
G=correlation_matrix(com,0);close(gcf);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 0, 400, 300]);imagesc(G(10:15,1:9));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:18]);yticks([1:1:18]);
% xticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});xtickangle(45)
% yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}','L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','aL2/3_{ex}','aL4_{ex}','aL2/3_{in}','aL4_{in}','med_{ex}','med_{in}','lat_{ex}','lat_{in}'});ylabel('Feature');
yticklabels({'PC1_{ex}','PC2_{ex}','PC3_{ex}','PC1_{in}','PC2_{in}','PC3_{in}'});xtickangle(45);set(gca,'FontSize',12)
xticklabels({'L2/3fr','L4fr','L2/3fr','L4fr','C alpha L23','C alpha L23','C X L23','C X L23','Pial depth'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',12)
%% Supplementary UMAP

close all;
%UMAP projections Panel C and D
%load embedding
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200508_embedding\reduced_data.mat')
plot_list = {'frac_vert','ang_inL23','ang_exL23','ang_inL4','ang_exL4','ang_inL5','ang_exL5','PCs'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%in vivo values non cyclic 

% %cyclic values
% plot_list ={'ORIpref','DIRpref'}
% plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'phasemap');
autoArrangeFigures;
%% 

figure(1);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(4);c=colorbar;caxis([0 0.6]);c.Ticks=[0:0.3:0.6];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(7);c=colorbar;caxis([0 85]);c.Ticks=[0:40:85];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(8);c=colorbar;caxis([0 140]);c.Ticks=[0:70:140];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(9);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(10);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(11);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(12);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(13);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(14);c=colorbar;caxis([0 1000]);c.Ticks=[0:500:1000];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(8);c=colorbar;caxis([0 140]);c.Ticks=[0:70:140];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(9);c=colorbar;caxis([0 160]);c.Ticks=[0:80:160];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(10);c=colorbar;caxis([100 320]);c.Ticks=[100:80:320];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(11);c=colorbar;caxis([0 130]);c.Ticks=[0:65:130];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(12);c=colorbar;caxis([50 350]);c.Ticks=[50:100:350];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(13);c=colorbar;caxis([0 500]);c.Ticks=[0:250:500];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(14);c=colorbar;caxis([200 500]);c.Ticks=[200:150:500];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(15);c=colorbar;caxis([0 500]);c.Ticks=[0:250:500];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(16);c=colorbar;caxis([200 500]);c.Ticks=[200:150:500];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(17);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(18);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(19);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(20);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(21);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(22);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;

% figure(15);c=colorbar;caxis([0 200]);c.Ticks=[0:100:200];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(16);c=colorbar;caxis([0 54]);c.Ticks=[0:27:54];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(17);c=colorbar;caxis([0 124]);c.Ticks=[0:62:124];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(18);c=colorbar;caxis([0 34]);c.Ticks=[0:17:34];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(19);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(20);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(21);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(22);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(23);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
% figure(24);c=colorbar;caxis([-2 1.7]);c.Ticks=[-2:1.5:1.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
close(figure(2));close(figure(3));close(figure(5));close(figure(6));
%% 
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp6\200511\'
savepdf_SW(fn,1);

%% In vivo parameters
close all;
 plot_list = {'OSIpref','DSIpref','ODIpref','ORIpref','DIRpref','Sigmapref','Capeakpref','sad','pci','noise','PCs','frac_vert'};
 plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
autoArrangeFigures;

figure(6);c=colorbar;caxis([2 4.3]);c.Ticks=[2:1:4.3];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(7);c=colorbar;caxis([1 6]);c.Ticks=[1:2:6];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(8);c=colorbar;caxis([0 0.16]);c.Ticks=[0:0.08:0.16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(9);c=colorbar;caxis([-4 1]);c.Ticks=[-4:2.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
figure(10);c=colorbar;caxis([-9.4 -2.4]);c.Ticks=[-9:3:-2.4];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;

%% Save

fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp6\200511\in vivo\'
savepdf_SW(fn,1);
%% Supplementary for Figure 5
close all;

%Rolling average plots
parameter_vector = abs(ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector);
ylabel('C_{x} L23 (µm)','Color','b');
ylim([5 60]);yticks([10:25:60]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

parameter_vector = abs(ang_exL23(:,3)-ang_exL23(:,1))*69;
rolling_avg_display(str,parameter_vector);
ylabel('C_{x} L23 (µm)','Color','r');
ylim([5 60]);yticks([10:25:60]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

parameter_vector = abs(ang_inL4(:,3)-ang_inL4(:,1))*69;
rolling_avg_display(str,parameter_vector);
ylabel('C_{x} L4 (µm)','Color','b');
ylim([5 70]);yticks([10:30:70]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

parameter_vector = abs(ang_exL4(:,3)-ang_exL4(:,1))*69;
rolling_avg_display(str,parameter_vector);
ylabel('C_{x} L4 (µm)','Color','r');
ylim([5 70]);yticks([10:30:70]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

%Rolling average plots L23alpha
parameter_vector =abs(a_L23in);
rolling_avg_display(str,parameter_vector,45,1)
ylabel('C_{\alpha} L23 (deg)','Color','r');
ylim([5 90]);yticks([10:20:90]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

parameter_vector = 90-abs(ang_exL4(:,5));
rolling_avg_display(str,parameter_vector)
ylabel('C_{\alpha} L4 (deg)','Color','r');
ylim([0 20]);yticks([0:10:20]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);

%% 


%dual barplots L4X
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_exL4(a,3)-ang_exL4(a,1))*69;
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('EX','Color','r');ylabel('C_{x} L4 (µm)','Color','k');set(gca,'FontSize',10); 
ylim([0 80]);ylabel('C_{x} L4 (µm)','Color','k')

%IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL4(a,3)-ang_inL4(a,1))*69;
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');set(gca,'FontSize',10); 
;ylabel('C_{x} L4 (µm)','Color','k');set(gca,'FontSize',10); 
ylim([0 80]);;
%% 

%dual barplots L23alpha
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=90-abs(ang_exL4(a,5));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45);
xlabel('EX','Color','r');
ylabel('C_{\alpha} L4 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 40]);
%IN
a=find(od_out_iviv(:,1)>0.25);  
par=90-abs(ang_inL4(a,5));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');
ylabel('C_{\alpha} L4 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 40]);

%dual barplots for L4 fraction in
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,1);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('L4fr','Color','r');set(gca,'FontSize',10);xtickangle(45)
%dual barplots for L4 fraction EX
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,2);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
;ylabel('L4fr','Color','b');set(gca,'FontSize',10);xtickangle(45)
%pia depth
a=find(od_out_iviv(:,1)>0.25);  
par=pia_input(a);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45)
%Correlation L4fr and CoMx/CoMalpha
a=find(od_out_iviv(:,1)>0.25); 
corr_plot(90-abs(ang_inL23(a,5)),L4fr(a,1),[],{'','',''});xlabel('L23 CoM\alpha (deg)','Color','b');
ylabel('L4fr','Color','r');set(gca,'FontSize',10);
%% 
a=find(od_out_iviv(:,1)>0.25);  
par=morph_sel(:,10);
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);xtickangle(45)

%% 
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp_rolling\Final_panels\'
savepdf_SW(fn,1);
%% Plot correlation for all pial depth and ang L4 L5 etc
%Pial depth vs CoMalpha L23 ex, 
corr_plot(90-abs((ang_exL23(:,5))),pia_input,[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('EX C_{\alpha} L2/3 (deg)','Color','r');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
%ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])

a=find(~isnan(ang_exL4(:,5))); 
corr_plot(90-abs((ang_exL4(a,5))),pia_input(a),[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('EX C_{\alpha} L4 (deg)','Color','r');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
%ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])

a=find(~isnan(ang_exL5(:,5))); 
corr_plot(90-abs((ang_exL5(a,5))),pia_input(a),[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('EX C_{\alpha} L5 (deg)','Color','r');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
%ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])

a=find(~isnan(ang_inL4(:,5))); 
corr_plot(90-abs((ang_inL4(a,5))),pia_input(a),[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('IN C_{\alpha} L4 (deg)','Color','b');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
%ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])

a=find(~isnan(ang_inL5(:,5))); 
corr_plot(90-abs((ang_inL5(a,5))),pia_input(a),[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('IN C_{\alpha} L5 (deg)','Color','b');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
%ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])

%% Try outs 

 data_in=ang_inL23;
for m=1:length(data_in)
              if data_in(m,10)==1 
                  temp(:,m)=(90-abs(data_in(m,5)))-180
              elseif data_in(m,10)==2 
                  temp(:,m)=abs((90-abs(data_in(m,5)))-180);
                   elseif data_in(m,10)==3 
                       temp(:,m)=(90-abs(data_in(m,5)))*-1
                        elseif data_in(m,10)==4 
                            temp(:,m)=(90-abs(data_in(m,5)))
              end
end
   
angle_h=temp';
%% 

close all
plot_list = {'frac_vert','ang_inL23','ang_exL23','span','ang_inL4','ang_exL4','PCs'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
autoArrangeFigures;



%% 

rolling_avg_display(str,angle_h,2)
%% 
rolling_avg_display(str,90-abs(ang_inL23(:,5)),2)
%% 

rolling_avg_display(str,(ang_inL23(:,3)-ang_inL23(:,1))*69,65,2)
%% 
a=find(od_out_iviv(:,1)>0.25);  
par=abs(angle_h(a));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
%% 
a=find(od_out_iviv(:,2)>0.25 & od_out_iviv(:,1)>0.25);  
par=(ang_inL23(a,3)-ang_inL23(a,1))*69
g1=find(od_out_iviv(a,5)>90 & od_out_iviv(a,5)<180) ;
g2=find(od_out_iviv(a,5)>275 & od_out_iviv(a,5)<360);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%xticklabels({'20-70°','100-150°'});

a=find(od_out_iviv(:,2)>0.25 & od_out_iviv(:,1)>0.25);  
par=(ang_exL23(a,3)-ang_exL23(a,1))*69
g1=find(od_out_iviv(a,5)>90 & od_out_iviv(a,5)<180) ;
g2=find(od_out_iviv(a,5)>270 & od_out_iviv(a,5)<360);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%xticklabels({'20-70°','100-150°'});

%% 
a=find(od_out_iviv(:,2)>0.25);
%corr_plot(od_out_iviv(a,4),90-abs(ang_inL23(a,5)),[],{'','',''});ylabel('C_{\alpha} (deg)','Color','b');
%xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0 180]);xticks([0:45:180]);
figure;polarscatter(90-abs(ang_inL23(a,5)),od_out_iviv(a,4))

rho=90-abs(ang_inL23(a,5))
theta=od_out_iviv(a,4);%or theta=linspace(0,180,numel(rho));
figure(1),clf(1)
polarscatter(deg2rad(theta),deg2rad(rho))
%% 
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
ylim([-50 100]);xlim([-100 100]);
%% 
a=find(od_out_iviv(:,2)>0.25);
fe=od_out_iviv(a,5);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'DIR'});
ylim([-50 100]);xlim([-100 100]);
%% 
pointsize=30
a=find(od_out_iviv(:,2)>0.25);
par1=(ang_exL23(a,3)-ang_exL23(a,1))*69
par2=(ang_inL23(a,3)-ang_inL23(a,1))*69
par3=od_out_iviv(a,2);
figure;
scatter(par1-par2,par3,'filled');
ylabel('DSI')
xlabel('EX-IN Cx')
%% 
%% 
rolling_avg_display(str,abs(angle_h),45,1)
ylabel('L2/3 alpha')
%% 
rolling_avg_display(str,angle_h,65,2)
ylabel('L2/3 alpha')
%% 
 data_in=ang_inL23;
for m=1:length(data_in)
              if data_in(m,10)==1 
                  temp(:,m)=(90-abs(data_in(m,5)))*-1
              elseif data_in(m,10)==2 
                  temp(:,m)=(90-abs(data_in(m,5)));
                   elseif data_in(m,10)==3 
                       temp(:,m)=(90-abs(data_in(m,5)))*-1
                        elseif data_in(m,10)==4 
                            temp(:,m)=(90-abs(data_in(m,5)));
              end
end
angle_ha=temp;
%% 

rolling_avg_display(str,angle_ha,65,2)
ylabel('L2/3 alpha')
%% 
par1=(ang_exL23(:,3)-ang_exL23(:,1))*69
par2=(ang_inL23(:,3)-ang_inL23(:,1))*69
rolling_avg_display(str,par1-par2,65,2)
%% 
%% 
rolling_avg_display(str,abs(scores(:,1)),45,1)
%% 
rolling_avg_display(str,180-abs(a_L23in_mod(:)),45,1)
%% 
a=find(od_out_iviv(:,1)>0.25);  
par=180-abs(a_L23in_mod(a));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});
%% 
[out_ang_inL23] = centroid_map(in_map(1:5,:,:),somac(:,1),somac(:,2),[1:147],2);
%% 
%% 
[a_L23in_mod]=sign_ang(out_ang_inL23);
%% 

corr_plot(180-abs(a_L23in_mod(a)),pia_input(a),[],{'','',''});ylabel('C_{\alpha} (deg)','Color','b');;xlabel('ORI')
%% 

%% Discussion 
a=find(od_out_iviv(:,1)>0.25);  
par=od_out_iviv(a,9)
g1=find(od_out_iviv(a,4)>15 & od_out_iviv(a,4)<65) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});

figure;histogram(par)

figure;scatter((ang_exL23(a,3)-ang_exL23(a,1))*69,(ang_inL23(a,3)-ang_inL23(a,1))*69)
%% 
a=find(od_out_iviv(:,1)>0.25);  
corr_plot(od_out_iviv(a,9),od_out_iviv(a,4),ang_inL23(a,8),{'','',''});ylabel('C_{\alpha} (deg)','Color','b');;xlabel('ORI')

%% 
 [od spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(str,1:147)
 %% 
 a=find(sftf_out_sel_iviv(:,2)>0.0);  
corr_plot(sftf_out_pref_iviv(a,1),ang_inL23(a,8),od_out_iviv(a,4),{'','',''});ylabel('C_{\alpha} (deg)','Color','b');;xlabel('ORI')
%% 
c_prf=vertcat(str(:).Ori)
diffb=abs(c_prf(bino_id,1)-c_prf(bino_id,2))

corr_plot(od_out_iviv(ipsi_id,4),ang_inL23(ipsi_id,8),ang_inL23(ipsi_id,8),{'','',''});
