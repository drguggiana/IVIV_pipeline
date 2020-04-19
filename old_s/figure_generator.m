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
ex_map_raw(:,:,i) = str(i).subpixel_raw_excMap;
in_map_raw(:,:,i) = str(i).subpixel_raw_inhMap;
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
% Setup A and Setup B
setups=[zeros(47,1);ones(100,1)];
% Slice orientation
slice_ori=[str(:).sliceOri];
%Fractions
frv=reshape([str(:).frac_vert],32,length(str))';
frh=reshape([str(:).frac_horz],32,length(str))';
L23fr=[sum(frv(:,3:5),2) sum(frv(:,19:21),2)];
L4fr=[sum(frv(:,6:7),2) sum(frv(:,22:23),2)];
L5fr=[sum(frv(:,8:10),2) sum(frv(:,24:26),2)];
%Differences of fractions
L5frt=L5fr;
L5frt(find(L5fr(:,1)==0),1)=NaN ;
L5frt(find(L5fr(:,2)==0),2)=NaN ;
diffL23fr=L23fr(:,1)-L23fr(:,2);
diffL4fr=L4fr(:,1)-L4fr(:,2);
diffL5fr=L5frt(:,1)-L5frt(:,2);
%span
span=reshape([str(:).span],6,length(str))';
%Principal compponent scores
scores=reshape([str(:).PCs],6,length(str))';
%subpixel soma
somac=reshape([str(:).subpixel_soma],2,length(str))';
%Centroid measurements
ang_exL23=reshape([str(:).ang_exL23],10,length(str))';
ang_exL4=reshape([str(:).ang_exL4],10,length(str))';
ang_exL5=reshape([str(:).ang_exL5],10,length(str))';
ang_inL23=reshape([str(:).ang_inL23],10,length(str))';
ang_inL4=reshape([str(:).ang_inL4],10,length(str))';
ang_inL5=reshape([str(:).ang_inL5],10,length(str))';
%Pial depth
pia_input=[str(:).pialD]';
%Orientations
oris=[0:45:315];
%visual responses properties
%1:OSI,2:DSI,3:ODI,4:ORI,5:DIR,6:Capeak,7:Tuning Width, 8:SF,9
%TF,10:SAD,11:Noise,12:PCI
od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
             ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci]]';
         
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
%EX vs IN horizontal centroid CoMx for L2/3, L4, L5, Panel C
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 500, 675, 200]);
subplot(1,3,1);
plot(ang_exL23(:,3)*69-ang_exL23(:,1)*69,ang_inL23(:,3)*69-ang_inL23(:,1)*69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];
ref.Color='k';box off;xlim([-150 150]);ylim([-150 150]);hold on;title('L23','FontWeight','normal');xticks([-150:75:150]);yticks([-150:75:150]);
subplot(1,3,2);
plot(ang_exL4(:,3)*69-ang_exL4(:,1)*69,ang_inL4(:,3)*69-ang_inL4(:,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];
ref.Color='k';box off;xlim([-150 150]);ylim([-150 150]);hold on;title('L4','FontWeight','normal');xticks([-150:75:200]);yticks([-150:75:150]);
subplot(1,3,3);
plot(ang_exL5(:,3)*69-ang_exL5(:,1)*69,ang_inL5(:,3)*69-ang_inL5(:,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');xlabel('EX','Color','r');ylabel('IN','Color','b');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;xlim([-300 300]);ylim([-300 300]);hold on;title('L5','FontWeight','normal');xticks([-300:150:300]);yticks([-300:150:300]);
%Pial depth vs CoMalpha L23 IN, Panel D
corr_plot(90-abs((ang_inL23(:,5))),pia_input,[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('CoM\alpha (deg)','Color','b');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
ylim([100 400]);yticks([100:150:400]);
%Pial depth vs PC1ex, Panel F
corr_plot(scores(:,1),pia_input,[],{'','',''});xlabel('PC1_{ex}','Color','r');
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);set(gca,'Ydir','reverse');ylim([100 400]);yticks([100:150:400]);
%CoMX vs PC13in, Panel G
corr_plot((ang_inL23(:,3)-ang_inL23(:,1))*69,scores(:,6),[],{'','',''});ylabel('PC3_{in}','Color','b');
xlabel('CoMx (µm)','Color','b');set(gca,'FontSize',10);
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure3\Final_panels\'
savepdf_SW(fn,1);
%% Figure 4 panels
close all;
%UMAP projections Panel C and D
%load embedding
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\embedding\reduced_data.mat');
%in vitro values
plot_list = {'frac_vert','ang_inL23','pialD'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%in vivo values
plot_list = {'OSIpref','DSIpref','ODIpref','ORIpref','DIRpref','Sigmapref','Capeakpref'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
autoArrangeFigures;

%Correlation Matrix Panel B
com=[];com=[L23fr(:,1) L4fr(:,1)  L23fr(:,2) L4fr(:,2) abs(ang_exL23(:,5)) abs(ang_inL23(:,5)) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) pia_input od_out_iviv(:,[1 2 3 4 5 6 7])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data OSI
a=find(od_out_iviv(:,1)>0.25)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) abs(ang_exL23(a,5)) abs(ang_inL23(a,5)) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
M=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data DSI
a=find(od_out_iviv(:,2)>0.25)
com=[];com=[L23fr(a,1) L4fr(a,1)  L23fr(a,2) L4fr(a,2) abs(ang_exL23(a,5)) abs(ang_inL23(a,5)) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
O=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
close(gcf);
%Assemble matrix
fG=[G([10 11 12],1:9) ; M([13],1:9) ; O([14],1:9) ; G([15 16],1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 300])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]);
xticklabels({'L2/3_{ex}','L4_{ex}','L2/3_{in}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});set(gca,'FontSize',10);
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%Selection of correlation from Panel B shown in Panel E
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),90-abs(ang_inL23(a,5)),[],{'','',''});ylabel('CoMx (µm)','Color','b');
xlabel('Orientation (deg)');set(gca,'FontSize',10);
%Selection of correlation from Panel B shown in Panel E
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),L4fr(a,1),[],{'','',''});ylabel('L4fr','Color','r');
xlabel('Orientation (deg)');set(gca,'FontSize',10);
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\Final_panels\'
savepdf_SW(fn,1);
%% Figure 5 panels
close all;
%Overview of centroid x and y with respect to soma for EX and IN and OSI colurcoded, Panel B
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
%Rolling average plots
parameter_vector = 90-abs(ang_inL23(:,5));
rolling_avg_display(str,parameter_vector)
ylabel('L23 CoM\alpha (µm)','Color','b');
ylim([10 70]);
set(gca,'FontSize',10);
%% 
%dual barplots
a=find(od_out_iviv(:,1)>0.25)  
par=90-abs(ang_inL23(a,5));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) 
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150) 
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2])
xticklabels({'20-90°','100-150°'})
xlabel('Orientation');
ylabel('CoM L23 (µm)','Color','b');  set(gca,'FontSize',12);
set(gca,'FontSize',12);