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
%ipsi, contra, bino, unres
ipsi_id=find([str(:).ipsi]==1);
contra_id=find([str(:).contra]==1);
bino_id=find([str(:).bino]==1);
unres_id=find([str(:).unres]==1);
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
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-100 100];
ref.YData=[-100 100];ylabel('IN C_{x} (µm)','Color','b')
ref.Color='k';box off;xlim([-100 100]);ylim([-100 100]);hold on;title('L23','FontWeight','normal');xticks([-100:50:100]);yticks([-100:50:100]);
subplot(1,3,2);
plot(ang_exL4(:,3)*69-ang_exL4(:,1)*69,ang_inL4(:,3)*69-ang_inL4(:,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');xlabel('EX C_{x} (µm)','Color','r');;ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];
ref.Color='k';box off;xlim([-150 150]);ylim([-150 150]);hold on;title('L4','FontWeight','normal');xticks([-150:75:200]);yticks([-150:75:150]);
subplot(1,3,3);
plot(ang_exL5(:,3)*69-ang_exL5(:,1)*69,ang_inL5(:,3)*69-ang_inL5(:,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;xlim([-300 300]);ylim([-300 300]);hold on;title('L5','FontWeight','normal');xticks([-300:150:300]);yticks([-300:150:300]);
%Pial depth vs CoMalpha L23 IN, Panel D
corr_plot(90-abs((ang_inL23(:,5))),pia_input,[],{'','',''});
ylabel('Pial depth (µm)','Color','k');xlabel('IN C_{\alpha} L2/3 (deg)','Color','b');set(gca,'FontSize',10);set(gca,'Ydir','reverse');
ylim([100 400]);yticks([100:150:400]);xlim([-5 100]);%xticks([-5:50:100])
%Pial depth vs PC1ex, Panel F
corr_plot(scores(:,1),pia_input,[],{'','',''});xlabel('PC1_{ex}','Color','r');
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);set(gca,'Ydir','reverse');ylim([100 400]);yticks([100:150:400]);xticks([-2:1:2]);
%CoMX vs PC13in, Panel G
corr_plot((ang_inL23(:,3)-ang_inL23(:,1))*69,scores(:,6),[],{'','',''});ylabel('PC3_{in}','Color','b');
xlabel('IN C_{x} L2/3 (µm)','Color','b');set(gca,'FontSize',10);ylim([-1.2 1.2]);yticks([-1.2:0.6:1.2]);xlim([-120 120]);xticks([-120:60:120])
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
plotting_embedding_str(reduced_data, str, plot_list, 0,1, 'parula');
%in vivo values non cyclic
plot_list = {'OSIpref','DSIpref','ODIpref','Sigmapref','Capeakpref'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%cyclic values
plot_list ={'ORIpref','DIRpref'}
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'phasemap');
autoArrangeFigures;

%Correlation Matrix Panel B
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) abs(ang_exL23(:,5)) abs(ang_inL23(:,5)) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) pia_input od_out_iviv(:,[1 2 3 4 5 6 7])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data OSI
a=find(od_out_iviv(:,1)>0.25)
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1) L4fr(a,2) abs(ang_exL23(a,5)) abs(ang_inL23(a,5)) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
M=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
close(gcf);
%Subsample data DSI
a=find(od_out_iviv(:,2)>0.25)
com=[];com=[L23fr(a,1) L23fr(a,2) L4fr(a,1)  L4fr(a,2) abs(ang_exL23(a,5)) abs(ang_inL23(a,5)) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) pia_input(a) od_out_iviv(a,[1 2 3 4 5 6 7])];
O=correlation_matrix(com,0);title('');xticks([1:1:17]);yticks([1:1:17]);
close(gcf);
%Assemble matrix
fG=[G([10 11 12],1:9) ; M([13],1:9) ; O([14],1:9) ; G([15 16],1:9)];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 200])
imagesc(fG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]);
xticklabels({'L2/3_{ex}','L2/3_{in}','L4_{ex}','L4_{in}','\alpha L2/3_{ex}','\alpha L2/3_{in}','COM L2/3_{ex}','COM L2/3_{in}','Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}'});set(gca,'FontSize',10);
c.Label.String = 'R';set(gca,'FontSize',12); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%Selection of correlation from Panel B shown in Panel E
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),90-abs(ang_inL23(a,5)),[],{'','',''});ylabel('C_{\alpha} (deg)','Color','b');
xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0 180]);xticks([0:45:180]);
%Selection of correlation from Panel B shown in Panel E
a=find(od_out_iviv(:,1)>0.25);
corr_plot(od_out_iviv(a,4),L4fr(a,1),[],{'','',''});ylabel('L4fr','Color','r');
xlabel('Orientation (deg)');set(gca,'FontSize',10);xlim([0 180]);xticks([0:45:180]);
%% 
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
close all;
%Overview of centroid x and y with respect to soma for EX and IN and OSI colurcoded, Panel B
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
%Rolling average plots L23alpha
parameter_vector = 90-abs(ang_inL23(:,5));
rolling_avg_display(str,parameter_vector)
ylabel('C_{\alpha} L23 (deg)','Color','b');
ylim([5 70]);set(gca,'FontSize',10);

%Rolling average plots
parameter_vector = abs(ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector);
ylabel('C_{\x} L23 (µm)','Color','b');
ylim([5 60]);set(gca,'FontSize',10);

%L4 angle alpha 
parameter_vector = 90-abs(ang_inL4(:,5));
rolling_avg_display(str,parameter_vector);
ylabel('L4 CoM\alpha (deg)','Color','b');


%dual barplots L23alpha
%EX
a=find(od_out_iviv(:,1)>0.25);  
par=90-abs(ang_exL23(a,5));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45);
xlabel('EX','Color','r');
ylabel('C_{\alpha} L23 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 100]);
%IN
a=find(od_out_iviv(:,1)>0.25);  
par=90-abs(ang_inL23(a,5));
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<70) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'20-70°','100-150°'});xtickangle(45)
xlabel('IN','Color','b');
ylabel('C_{\alpha} L23 (deg)','Color','k'); ;set(gca,'FontSize',10);
ylim([0 100]);
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
%plot peak nr of sholl crossings BASAL vs TW, PANEL E
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 30]);
ylabel('Tuning Width','Color','k');xlabel('Peak Nr. sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%Plot correlation matrix, PANEL F
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[morph_parameters(:,2) nanmax(df,[],2) dis_s'  morph_parameters(:,4)  max_s'  sum_densap' ...
    morph_parameters(:,12) nanmax(db,[],2) dis_s_ba'  morph_parameters(:,14) max_s_ba'  sum_densba'  od_out_iviv(:,[1 2 3 4 5 7 8]) pia_input]
G=correlation_matrix(com,0);close(gcf)
Gf=G(13:end,1:12)
Gf(6,1)=0;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(Gf);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:8]);
 
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Density','Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Density'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','ORI','DIR','TW','Ca_{peak}','Pial depth'});ytickangle(45);set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
% Plot the input vs morphology density,panel G
for i=1:length(str)
%whole map
tmp1=ex_map(3:end,:,i);
tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:));sum(tmp2(:))];
end
corr_plot(max_s(morph_cells_id)',tot_input(morph_cells_id,1),[],{'','',''});xlabel('Peak Nr. sholl crossings','Color','k');
ylabel('Total input','Color','r');
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure6\Final_panels\'
savepdf_SW(fn,1);
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
b=od_out(:,1)>0.25;
subplot(1,2,2);
binRange = 0:45:179;
hcx = histcounts(od_out(b,4),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,4),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Orientation (deg)','FontSize',12);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
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
b=od_out(:,2)>0.25;
subplot(1,2,2);
binRange = 0:45:315;
hcx = histcounts(od_out(b,5),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,5),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Direction (deg)','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]); 
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
    length(ipsi_id)/resp_iviv length(bino_id)/resp_iviv],'stacked');box off;ylabel('Fraction cells');
 b(1, 1).FaceColor='b';b(1, 2).FaceColor='r';b(1, 3).FaceColor='w';
 xticklabels({'all','iviv'});legend('contra','ipsi','bino');legend boxoff;set(gca,'FontSize',12);
%TW distribution 
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 250, 200]);set(gcf,'color','w');
binRange = 0:20:90;
hcx = histcounts(od_out(find(~isnan(od_out(:,7))),7),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,7))),7),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Tuning width (deg)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.6]);yticks([0:0.2:0.6]);
%Spon activity
fig23 = figure;
set(fig23, 'Name', 'Binocular protocol');set(fig23, 'Position', [1000, 800, 600, 200]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0.0013:0.025:0.2;
hcx = histcounts(spon_out(find(~isnan(spon_out(:,1))),1),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,10))),10),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('SAD (Hz)','FontSize',12);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];

ylim([0 0.8]);yticks([0:0.4:0.8]);
subplot(1,2,2);
binRange = -7:1:2;
hcx = histcounts(log(abs(od_out_iviv(:,12))),[binRange Inf],'Normalization','probability');
hcy = histcounts(log(abs(spon_out(:,2))),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('PCI','FontSize',12);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.8]);yticks([0:0.4:0.8]);
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp3\Final_panels\'
savepdf_SW(fn,1);