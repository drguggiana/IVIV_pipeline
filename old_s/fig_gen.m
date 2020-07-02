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

frh=reshape([str(:).frac_horz],32,length(str))';
frv=reshape([str(:).frac_vert],32,length(str))';
L1fr=[sum(frv(:,17:18),2)];
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
ang_exL23=reshape([str(:).ang_exL23],10,length(str))';
ang_exL4=reshape([str(:).ang_exL4],10,length(str))';
ang_exL5=reshape([str(:).ang_exL5],10,length(str))';
ang_inL23=reshape([str(:).ang_inL23],10,length(str))';
ang_inL4=reshape([str(:).ang_inL4],10,length(str))';
ang_inL5=reshape([str(:).ang_inL5],10,length(str))';

ang_maxL23ex=reshape([str(:).max_ex],length(str),9);

ang_maxL23in=reshape([str(:).max_in],length(str),9);

ang_exL23nw=reshape([str(:).ang_exL23_nonweighted],10,length(str))';
ang_inL23nw=reshape([str(:).ang_inL23_nonweighted],10,length(str))';
%Pial depth
pia_input=[str(:).pialD]';
%Orientations
oris=[0:45:315];
od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
             ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci];[str(:).error_pref];[str(:).r2_pref];[str(:).tun_pref]]';
%ipsi, contra, bino, unres, resp
 ipsi_id=find([str(:).ipsi]==1);
 contra_id=find([str(:).contra]==1);
 bino_id=find([str(:).bino]==1);
 unres_id=find([str(:).unres]==1);
 resp_id=find([str(:).resp]==1);
 %% Plot data with fits and calculate R2
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 200, 1600, 800]);
for i=1:length(resp_id)
subplot(6,10,i)
hold on;plot(oris,str(resp_id(i)).TCpref,'--o');
y=str(resp_id(i)).TCpref;
ym=mean(y)
;xticks([0:90:360]);
hold on;
if str(resp_id(i)).ipsi==1 
 plot(str(resp_id(i)).fit_resp(361:end));
 yp=(str(resp_id(i)).fit_resp(361:end));
 yp_s=yp(oris+1);
 
elseif str(resp_id(i)).contra==1 
     plot(str(resp_id(i)).fit_resp(1:360));
     yp=(str(resp_id(i)).fit_resp(1:360));
 yp_s=yp(oris+1);
else str(resp_id(i)).bino==1 
    if str(resp_id(i)).ODIpref>0 
    hold on;plot(str(resp_id(i)).fit_resp(1:360))
        yp=(str(resp_id(i)).fit_resp(1:360));
 yp_s=yp(oris+1);
    else str(resp_id(i)).ODIpref<0 
        ;hold on;plot(str(resp_id(i)).fit_resp(361:end))
            yp=(str(resp_id(i)).fit_resp(361:end));
 yp_s=yp(oris+1);
    end
end
r2=1-(sum(sqrt((y-yp_s').^2))/sum(sqrt((y-ym).^2)))
r_square(i)=r2
title([num2str(str(resp_id(i)).OSIpref) ' / ' num2str(r_square(i))]);
y=[]
ym=[];
yp=[];
yp_s=[];
r2=[];
end
%r2 = 1-(sum(sqrt((y-yp).^2))/sum(sqrt((y-ym).^2)));
r_sq=ones(147,1)*NaN;
r_sq(resp_id)=r_square;
%% Add r2 square to str
for i=1:length(r_sq)
    str(i).r_sq=r_sq(i);
end
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

%EX vs IN horizontal centroid CoMx for L2/3, L4, L5, Panel C
g=find(ang_inL4(:,3)*69-ang_inL4(:,1)*69<135)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 675]);
subplot(3,2,1);
plot(ang_exL23(:,3)*69-ang_exL23(:,1)*69,ang_inL23(:,3)*69-ang_inL23(:,1)*69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-100 100];
ref.YData=[-100 100];
ref.Color='k';box off;xlim([-115 100]);ylim([-115 100]);hold on;title('L23','FontWeight','normal');xticks([-100:50:100]);yticks([-100:50:100]);
subplot(3,2,3);
plot(ang_exL4(g,3)*69-ang_exL4(g,1)*69,ang_inL4(g,3)*69-ang_inL4(g,1)*69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];ylabel('IN C_{x} (µm)','Color','b')
ref.Color='k';box off;xlim([-175 150]);ylim([-175 150]);hold on;title('L4','FontWeight','normal');xticks([-150:75:200]);yticks([-150:75:150]);
subplot(3,2,5);
plot(ang_exL5(:,3)*69-ang_exL5(:,1)*69,ang_inL5(:,3)*69-ang_inL5(:,1)*69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;xlim([-335 300]);ylim([-335 300]);hold on;title('L5','FontWeight','normal');xticks([-300:150:300]);yticks([-300:150:300]);
xlabel('EX C_{x} (µm)','Color','r')
ref.XData=[-300 300];
ref.YData=[-300 300];
%EX vs IN verticalcentroid CoMY for L2/3, L4, L5, Panel C
subplot(3,2,2);
plot(ang_exL23(:,4)*-69-ang_exL23(:,2)*-69,ang_inL23(:,4)*-69-ang_inL23(:,2)*-69,'.','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];;xlim([-175 150]);ylim([-175 150]);xticks([-150:100:150]);yticks([-150:100:150]);
ref.Color='k';box off;;hold on;title('L23','FontWeight','normal');
subplot(3,2,4);
plot(((ang_exL4(g,4)*-69))-ang_exL4(g,2)*-69,((ang_inL4(g,4)*-69))-ang_inL4(g,2)*-69,'.','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',10);
set(gcf,'color','w');;xlim([-365 -25]);ylim([-365 -25]);
;ref= refline(1,0);
set(gca,'FontSize',10);ref.LineStyle='--';
ref.XData=[-340 -25];
ref.YData=[-340 -25];
ylabel('IN C_{y} (µm)','Color','b')
ref.Color='k';box off;hold on;title('L4','FontWeight','normal');xticks([-325:100:-25]);yticks([-325:100:-25]);
subplot(3,2,6);
plot(((ang_exL5(:,4)*-69))-ang_exL5(:,2)*-69,((ang_inL5(:,4)*-69))-ang_inL5(:,2)*-69,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--';ref.XData=[-500 -50];
ref.YData=[-500 -50];xlim([-525 -50]);ylim([-525 -50]);xticks([-500:225:-50]);yticks([-500:225:-50]);
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
%% Principal component corr matrix
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
corr_plot(scores(:,5),abs((ang_inL23(:,4)-ang_inL23(:,2))*-69),[],{'','',''});xlabel('PC2_{in}','Color','b');xlim([-2.2 2]);xticks([-2:1:2]);ylim([-10 150]);yticks([0:50:150])
%PC1ex with 
corr_plot(scores(:,1),pia_input,[],{'','',''});xlabel('PC1_{ex}','Color','r');
ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10);set(gca,'Ydir','reverse');ylim([100 420]);yticks([100:150:400]);xticks([-2:1:2]);xlim([-2.2 2])
%CoMX vs PC13in, Panel G
corr_plot(scores(:,6),(ang_inL23(:,3)-ang_inL23(:,1))*69,[],{'','',''});xlabel('PC3_{in}','Color','b');
ylabel('IN C_{x} L2/3 (µm)','Color','b');set(gca,'FontSize',10);xlim([-1.4 1.2]);xticks([-1.2:0.6:1.2]);ylim([-140 120]);yticks([-120:60:120])

%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure3\Final_panels\'
savepdf_SW(fn,1);


%% Figure 4 panels
%Example cells
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
%% UMAP
close all;
%UMAP projections Panel C and D
%load embedding
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200523_embedding\reduced_data.mat');
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
%% Rearrange umap plots and scale
 %figure(1);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];
 figure(2);c=colorbar;caxis([0 0.7]);c.Ticks=[0:0.35:0.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(2), 'Position', [200, 0, 200, 200])
 figure(4);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(3), 'Position', [200, 0, 200, 200])
 figure(7);c=colorbar;caxis([0 120]);c.Ticks=[0:60:120];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(5), 'Position', [200, 0, 200, 200])
 figure(8);c=colorbar;caxis([0 140]);c.Ticks=[0:70:140];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(6), 'Position', [200, 0, 200, 200])
 figure(9);c=colorbar;caxis([120 360]);c.Ticks=[120:120:360];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(7), 'Position', [200, 0, 200, 200])
% figure(8);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;
 figure(16);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(14), 'Position', [200, 0, 200, 200])
figure(17);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(15), 'Position', [200, 0, 200, 200])
 figure(18);c=colorbar;caxis([-1 1]);c.Ticks=[-1:1:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(16), 'Position', [200, 0, 200, 200])
 figure(19);c=colorbar;caxis([2 4]);c.Ticks=[2:5];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(17), 'Position', [200, 0, 200, 200])
 figure(21);c=colorbar;caxis([0 180]);c.Ticks=[0:45:180];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(19), 'Position', [200, 0, 200, 200])
 figure(22);c=colorbar;caxis([0 360]);c.Ticks=[0:90:360];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(20), 'Position', [200, 0, 200, 200])
%  close(figure(1));close(figure(4));close(figure(9));close(figure(10));close(figure(11));close(figure(12));close(figure(13));
%  close(figure(8));close(figure(17));close(figure(18));
%% Save umap figures
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200523_embedding\'
savepdf_SW(fn,1);
%% Creating correlation matrix
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) L5fr(:,1) L5fr(:,2) abs(ang_exL23(:,3)-ang_exL23(:,1))...
    abs(ang_inL23(:,3)-ang_inL23(:,1)) abs(ang_exL23(:,4)-ang_exL23(:,2))...
   abs(ang_inL23(:,4)-ang_inL23(:,2))  pia_input od_out_iviv(:,[1 2 3 6])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%% Subsample sigma based on quality of the fit
a=find(od_out_iviv(:,7)>0.3);
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) abs(ang_exL23(a,4)-ang_exL23(a,2))...
   abs(ang_inL23(a,4)-ang_inL23(a,2))  pia_input(a) od_out_iviv(a,[7])]; 
G2=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%% Circular correlation for ORI
a=find(od_out_iviv(:,1)>0.25); 
par_c=[];
 par_c=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1)  L5fr(a,2) abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) abs(ang_exL23(a,4)-ang_exL23(a,2))...
   abs(ang_inL23(a,4)-ang_inL23(a,2)) pia_input(a)]
for i=1:11
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(a,4)), par_c(:,i))
end
% Circular correlation for DRI
a=find(od_out_iviv(:,2)>0.25); 
par_c=[];
 par_c=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1)  L5fr(a,2)  abs(ang_exL23(a,3)-ang_exL23(a,1))...
    abs(ang_inL23(a,3)-ang_inL23(a,1)) abs(ang_exL23(a,4)-ang_exL23(a,2))...
   abs(ang_inL23(a,4)-ang_inL23(a,2)) pia_input(a)]
for i=1:11
    [rho2(i) pval2(i)] = circ_corrcl(deg2rad(od_out_iviv(a,5)), par_c(:,i))
end
%% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% Combine matrixes and Plot combined correlatioan matrix
fG=[];fG2=[];
fG=[G([12 13 14 15],1:11)];
fG2=[G2(12,1:11)]
fG=[fG; fG2;c_cor_e];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 200])
imagesc(fG);c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:12]);yticks([1:1:12])
yticklabels({'gOSI','gDSI','ODI','Ca_{peak}','TW*','ORI*','DIR*'});set(gca,'FontSize',10);
xticklabels({'L2/3','L2/3','L4','L4','L5','L5','C L2/3 ','C L2/3','C L2/3','C L2/3', 'Pial depth'});xtickangle(45);set(gca,'FontSize',10);
%% %% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\Final_panels\'
savepdf_SW(fn,1);

%% Figure 5 panels
%Overview of centroid x and y with respect to soma for EX and IN and OSI colurcoded, Panel B
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
%% Rolling average and barplots
%define two barplots windows
s1a=10;s1b=60;s2a=100;s2b=150
%% Rolling average  Cy
parameter_vector = abs(ang_exL23(:,4)*69-ang_exL23(:,2)*69)
parameter_vector2 = abs(ang_inL23(:,4)*69-ang_inL23(:,2)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L2/3 (µm)','Color','k');
ylim([5 70]);yticks([10:30:70]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10);
%% %% Barplots Cy  EX and IN, L4 and L5 fractions
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_exL23(a,4)*69-ang_exL23(a,2)*69)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylim([0 80]);yticks([0:40:80]);ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10); 
%IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL23(a,4)*69-ang_inL23(a,2)*69)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylim([0 80]);yticks([0:40:80]);ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);
%L4 EX
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,1)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','k'); ;set(gca,'FontSize',10); ylim([0 0.6]);yticks([0:0.3:0.6]);
%L4 IN
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,2)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','k'); ;set(gca,'FontSize',10); ylim([0 0.6]);yticks([0:0.3:0.6])
%L5 EX
a=find(od_out_iviv(:,1)>0.25);  
par=L5fr(a,1)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,2);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L5fr','Color','k'); ;set(gca,'FontSize',10); ylim([0 0.4]);yticks([0:0.2:0.4]);
%% Rolling average  Cx
parameter_vector = abs(ang_exL23(:,3)*69-ang_exL23(:,1)*69)
parameter_vector2 = abs(ang_inL23(:,3)*69-ang_inL23(:,1)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cx L2/3 (µm)','Color','k');
ylim([5 70]);yticks([10:30:70]);xlim([0 180]);xticks([0:45:180]);set(gca,'FontSize',10); 
%% Barplots Cx EX and IN
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);  
par=abs(ang_exL23(a,3)*69-ang_exL23(a,1)*69)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylim([0 110]);yticks([0:55:110]);ylabel('Cx L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10); 
%IN
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);  
par=abs(ang_inL23(a,3)*69-ang_inL23(a,1)*69)
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);ylim([0 110]);yticks([0:55:110]);ylabel('Cx L2/3 (µm)','Color','k') ;set(gca,'FontSize',10);
%% Correlation comparison across layers
g1=[];g2=[];
g1=find(od_out_iviv(:,1)>0.25 & ang_inL4(:,3)*69-ang_inL4(:,1)*69<135 & od_out_iviv(:,4)>s1a & od_out_iviv(:,4)<s1b);
g2=find(od_out_iviv(:,1)>0.25 & ang_inL4(:,3)*69-ang_inL4(:,1)*69<135 & od_out_iviv(:,4)>s2a & od_out_iviv(:,4)<s2b);
par1=[];par2=[];
par1=ang_exL23(g1,3)*69-ang_exL23(g1,1)*69
par2=ang_exL4(g1,3)*69-ang_exL4(g1,1)*69
[R1,P1,RLO1,RUP1]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_inL23(g1,3)*69-ang_inL23(g1,1)*69
par2=ang_inL4(g1,3)*69-ang_inL4(g1,1)*69
[R4,P4,RLO4,RUP4]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_exL23(g2,3)*69-ang_exL23(g2,1)*69
par2=ang_exL4(g2,3)*69-ang_exL4(g2,1)*69
[R7,P7,RLO7,RUP7]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_inL23(g2,3)*69-ang_inL23(g2,1)*69
par2=ang_inL4(g2,3)*69-ang_inL4(g2,1)*69
[R10,P10,RLO10,RUP10]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05);
%Figure
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 200, 300]);
scatter([R1(2); R4(2)],[1 3],'m','filled')
hold on;scatter([R7(2); R10(2)],[2 4],'c','filled')
ylim([0 5]);yticks([1:1:5]);xlabel('Correlation');
hold on;plot([RLO1(2) RUP1(2)],[1 1],'-m');
hold on;plot([RLO4(2) RUP4(2)],[3 3],'-m');
hold on;plot([RLO7(2) RUP7(2)],[2 2],'-c');
hold on;plot([RLO10(2) RUP10(2)],[4 4],'-c');
yticklabels({'EX','EX','IN','IN'});title('aCxL2/3 to aCxL4');
legend('10-60°','100-150');legend boxoff;
set(gca,'FontSize',10);
%% Rolling average  vector length
parameter_vector = ang_exL23(:,8)*69
parameter_vector2 = ang_inL23(:,8)*69
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cd L2/3 (µm)','Color','k');
   ylim([10 90]);yticks([10:40:90]);
 xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10)
%% Save figures as PDF in a folder specified by fn IF DESIRED
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure5\Final_panels\'
savepdf_SW(fn,1);

%% Figure 6
%% Perform sholl anaylsis, this should be incoorporated into the structure
%Get all morphtraces
close all;
for i=1:length(str)
    if ~isempty(str(i).morph)==1 
    zz{:,i}=str(i).morphtraces;
    else
        zz{:,i}=NaN;
    end
end
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147,1);
%% 
close all;
%Plot representative examples, cell 143 and 122 for panel A, 
%NOTE tuning curve: individual traces not in the str only the average atm
 id_m=143;
 plot_morphologies(str,id_m,1,1,1);
id_m2=122;
 plot_morphologies(str,id_m2,1,1,1);

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
%% %subsample cells with morph and fct. in vivo
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1 & str(i).resp==1 & ~isnan(str(i).OSIpref)==1
        m_res(i)=1;     
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=find(m_res==1);
%R2 criterion
morph_res_sub2=morph_res_sub;
morph_res_sub2(find(r_sq(morph_res_sub2)<0.3))=[];
%Plot correlation Wuning width
corr_plot(morph_parameters(morph_res_sub2,4),od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 25]);
ylabel('Tuning Width','Color','k');xlabel('Nr of branch points','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
corr_plot(max_s(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 15]);
ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%plot peak nr of sholl crossings BASAL vs TW, PANEL E
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 30]);
ylabel('Tuning Width','Color','k');xlabel('Peak Nr. sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%% Plot correlation matrix, PANEL F
df=[morph_parameters(morph_res_sub,9) morph_parameters(morph_res_sub,10)];
db=[morph_parameters(morph_res_sub,19) morph_parameters(morph_res_sub,20)];
com=[];com=[morph_parameters(morph_res_sub,2) nanmax(df,[],2) dis_s(morph_res_sub)'  morph_parameters(morph_res_sub,4)  max_s(morph_res_sub)' ...
    morph_parameters(morph_res_sub,12) nanmax(db,[],2) dis_s_ba(morph_res_sub)'  morph_parameters(morph_res_sub,14) max_s_ba(morph_res_sub)'  od_out_iviv(morph_res_sub,[1 2 3 6]) pia_input(morph_res_sub)]
G1=correlation_matrix(com,0);
%% subsample TW
df=[morph_parameters(morph_res_sub2,9) morph_parameters(morph_res_sub2,10)];
db=[morph_parameters(morph_res_sub2,19) morph_parameters(morph_res_sub2,20)];
com=[];com=[morph_parameters(morph_res_sub2,2) nanmax(df,[],2) dis_s(morph_res_sub2)'  morph_parameters(morph_res_sub2,4)  max_s(morph_res_sub2)' ...
    morph_parameters(morph_res_sub2,12) nanmax(db,[],2) dis_s_ba(morph_res_sub2)'  morph_parameters(morph_res_sub2,14) max_s_ba(morph_res_sub2)'  od_out_iviv(morph_res_sub2,[7])]
G2=[];
G2=correlation_matrix(com,0);
%% Circular correlation for ORI
a=[];
a=find(od_out_iviv(morph_res_sub,1)>0.25); 
par_c=[];rho1=[];pval1=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),4)), par_c(:,i))
end
% Circular correlation for DRI
a=[]
a=find(od_out_iviv(morph_res_sub,2)>0.25); 
par_c=[];
rho2=[];pval2=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho2(i) pval2(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),5)), par_c(:,i))
end
%% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% 
Gf=G1(11:14,1:10)
tG=[Gf;G2(11:end,1:10);c_cor_e;G1(15,1:10)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(tG);c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:8]) 
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'OSI','DSI','ODI','Ca_{peak}','TW*','ORI*','DIR*','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);








































%% Supplementary Figures
%% Cell plotter, morphology, in vivo, maps , SUPPLEMENTRAY 1
close all;
for i=1:147;
figure;
iviv_plotter(str,i)
end
%% Save individual morpho, in vivo, maps panels
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp1\Block3\'
savepdf_SW(fn,0);
%% Example for with and without TTX
pathName='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\TTX\161116\SW0005_inhibition_before\map04'
load_raw_map(pathName,str,68)
%% PCA supplementary: PCs eigenvector map plots, varaince expalined PCcom and C3x/in correlation matrix
close all;
map_align_PCA(str);
%% Correlation matrix Pial depth Cx Cy layers for supplement + PCs and real data corr matrix
com=[];com=[ang_exL23(:,3)-ang_exL23(:,1)...
    ang_inL23(:,3)-ang_inL23(:,1) ang_exL23(:,4)-ang_exL23(:,2) ang_inL23(:,4)-ang_inL23(:,2)... 
    ang_exL4(:,3)-ang_exL4(:,1) ang_inL4(:,3)-ang_inL4(:,1) ang_exL4(:,4)-ang_exL23(:,2) ang_inL4(:,4)-ang_inL4(:,2)... 
    ang_exL5(:,3)-ang_exL5(:,1) ang_inL5(:,3)-ang_inL5(:,1) ang_exL5(:,4)-ang_exL5(:,2) ang_inL5(:,4)-ang_inL5(:,2)...
    pia_input]; 
O=correlation_matrix(com,0);title('');

xticks([1:1:13]);yticks([1:1:13]);caxis([-1 1]);c=colorbar;c.Ticks=[-1:0.5:1];
xticklabels({'C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','PiaD'});set(gca,'FontSize',10);xlabel('');set(gca,'FontSize',12)
yticklabels({'C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','C_{x}','C_{x}','C_{y}','C_{y}','PiaD'});set(gca,'FontSize',10);ylabel('');set(gca,'FontSize',12)
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


%% in vivo in vitro comparison
% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 [od_out sftf_out sftf_sel sftf_pref spon_out pia_all delta_Ca fit_Ca] = concat_invivo(L23_PC);
%% calculate R2 of fit
yp=fit_Ca(:,oris+1);
for i=1:length(fit_Ca)
    y=delta_Ca(i,:)
    ym=mean(y);
    yp=fit_Ca(i,oris+1);
r2=1-(sum(sqrt((y-yp).^2))/sum(sqrt((y-ym).^2)))
r2_all(i)=r2;
end
 figure;set(gcf,'color','w');scatter(od_out(:,1),r2_all);box off;
 xlabel('gOSI');ylabel('R2 of fit');
%% OSI and ORI
close all;
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [200, 200, 600, 200]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(od_out(find(~isnan(od_out(:,1))),1),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,1))),1),[binRange Inf],'Normalization','probability');
b1=bar([0.1:0.2:0.9] ,[hcx(1:5);hcy(1:5)]');box off;xlabel('OSI','FontSize',12);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);
xlim([-0.1 1]);xticks([0.1:0.2:0.9])
a=od_out_iviv(:,1)>0.25;
b=od_out(:,1)>0.25;set(gca,'FontSize',10);xtickangle(45);
subplot(1,2,2);
binRange = [0 45 90 135 180 225]-22.5;
hcx = histcounts(od_out(b,4),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,4),[binRange Inf],'Normalization','probability');
hcy_m=[hcy(1)+hcy(5) hcy(2:4)];
hcx_m=[hcx(1)+hcx(5) hcx(2:4)];
b1=bar([0 45 90 135],[hcx_m;hcy_m]');box off;xlabel('Orientation (°)','FontSize',12);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);
%% DSI and DIR
fig21 = figure;
set(fig21, 'Name', 'Binocular protocol');set(fig21, 'Position', [800, 200, 600, 200]);set(gcf,'color','w');
subplot(1,2,1);
binRange = 0:0.2:1;
hcx = histcounts(od_out(find(~isnan(od_out(:,2))),2),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,2))),2),[binRange Inf],'Normalization','probability');
b1=bar([0.1:0.2:0.9],[hcx(1:5);hcy(1:5)]');box off;xlabel('DSI','FontSize',12);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.4]);yticks([0:0.2:0.4]);xlim([-0.1 1]);xticks([0.1:0.2:0.9])
a=od_out_iviv(:,2)>0.25;
b=od_out(:,2)>0.25;set(gca,'FontSize',10);
subplot(1,2,2);
binRange = [0:45:405]-22.5;
hcx = histcounts(od_out(b,5),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(a,5),[binRange Inf],'Normalization','probability');
hcy_m=[hcy(1)+hcy(9) hcy(2:8)];
hcx_m=[hcx(1)+hcx(9) hcx(2:8)];
b1=bar([0:45:315],[hcx_m;hcy_m]');box off;xlabel('Direction (deg)','FontSize',10);
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
xtickangle(45)
ylim([0 0.4]);yticks([0:0.2:0.4]);set(gca,'FontSize',10);
%% ipsi contra bino
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
 
  fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 300, 200]);set(gcf,'color','w');
h=histogram(od_out(find(~isnan(od_out(:,3))),3),7,'Normalization','probability');h.FaceColor=[0 0 0];
hold on;h=histogram(od_out_iviv(find(~isnan(od_out_iviv(:,3))),3),7,'Normalization','probability');
h.FaceColor=[0 0.7 0.6];box off;
ylim([0 0.3]);yticks([0:0.15:0.3]);xlabel('ODI');ylabel('Fraction of cells');set(gca,'FontSize',10);
%% TW
c=r2_all>0.3;set(gca,'FontSize',10);
d=r_sq>0.3;set(gca,'FontSize',10);
%TW distribution 
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 250, 200]);set(gcf,'color','w');
binRange = 0:20:90;
hcx = histcounts(od_out(c,7),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(d,7),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Tuning width (deg)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.6]);yticks([0:0.2:0.6]);set(gca,'FontSize',10);
%% Ca peak
fig23=figure;set(gcf,'color','w');set(fig23, 'Position', [800, 600, 250, 200]);set(gcf,'color','w');
binRange = 0:50:400;
hcx = histcounts(od_out(find(~isnan(od_out(:,6))),6),[binRange Inf],'Normalization','probability');
hcy = histcounts(od_out_iviv(find(~isnan(od_out_iviv(:,6))),6),[binRange Inf],'Normalization','probability');
b1=bar(binRange,[hcx;hcy]');box off;xlabel('Ca_{peak} (?R/R0)','FontSize',10);ylabel('Fraction of cells');
b1(1).FaceColor=[0 0 0];
b1(2).FaceColor=[0 0.7 0.6];
ylim([0 0.7]);yticks([0:0.35:0.7]);set(gca,'FontSize',10);xtickangle(45);
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

%% 
% tm1=length(find(sftf_out(:,1)==0.02 & sftf_out(:,2)==1))/sum(~isnan(sftf_out(:,1)));
% tm2=length(find(sftf_out(:,1)==0.02 & sftf_out(:,2)==2))/sum(~isnan(sftf_out(:,1)));
% tm3=length(find(sftf_out(:,1)==0.02 & sftf_out(:,2)==4))/sum(~isnan(sftf_out(:,1)));
% tm4=length(find(sftf_out(:,1)==0.08 & sftf_out(:,2)==1))/sum(~isnan(sftf_out(:,1)));
% tm5=length(find(sftf_out(:,1)==0.08 & sftf_out(:,2)==2))/sum(~isnan(sftf_out(:,1)));
% tm6=length(find(sftf_out(:,1)==0.08 & sftf_out(:,2)==4))/sum(~isnan(sftf_out(:,1)));
% tm7=length(find(sftf_out(:,1)==0.16 & sftf_out(:,2)==1))/sum(~isnan(sftf_out(:,1)));
% tm8=length(find(sftf_out(:,1)==0.16 & sftf_out(:,2)==2))/sum(~isnan(sftf_out(:,1)));
% tm9=length(find(sftf_out(:,1)==0.16 & sftf_out(:,2)==4))/sum(~isnan(sftf_out(:,1)));
%sftfall=[tm1 tm4 tm7;tm2 tm5 tm8;tm3 tm6 tm9];
tm1=length(find(od_out_iviv(:,8)==0.02 & od_out_iviv(:,9)==1))/sum(~isnan( od_out_iviv(:,9)));
tm2=length(find(od_out_iviv(:,8)==0.02 &  od_out_iviv(:,9)==2))/sum(~isnan( od_out_iviv(:,9)));
tm3=length(find(od_out_iviv(:,8)==0.02 &  od_out_iviv(:,9)==4))/sum(~isnan( od_out_iviv(:,9)));
tm4=length(find(od_out_iviv(:,8)==0.08 &  od_out_iviv(:,9)==1))/sum(~isnan( od_out_iviv(:,9)));
tm5=length(find(od_out_iviv(:,8)==0.08 &  od_out_iviv(:,9)==2))/sum(~isnan( od_out_iviv(:,9)));
tm6=length(find(od_out_iviv(:,8)==0.08 &  od_out_iviv(:,9)==4))/sum(~isnan( od_out_iviv(:,9)));
tm7=length(find(od_out_iviv(:,8)==0.16 &  od_out_iviv(:,9)==1))/sum(~isnan( od_out_iviv(:,9)));
tm8=length(find(od_out_iviv(:,8)==0.16 &  od_out_iviv(:,9)==2))/sum(~isnan( od_out_iviv(:,9)));
tm9=length(find(od_out_iviv(:,8)==0.16 &  od_out_iviv(:,9)==4))/sum(~isnan( od_out_iviv(:,9)));
sftfalli=[tm1 tm4 tm7;tm2 tm5 tm8;tm3 tm6 tm9];
%% 

fig23 = figure;
set(fig23, 'Name', 'SFTF');set(fig23, 'Position', [200, 800, 400, 200]);
set(gcf,'color','w');
subplot(1,2,1);imagesc(sftfall);c=colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);caxis([0 0.25]);
set(gca,'YTickLabel',{'1' '2','4'});set(gca,'Ydir','reverse');ylabel('TF');
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;xtickangle(45);
hold on; title(['all' ' n=' num2str(sum(~isnan(sftf_out(:,1))))]);caxis([0 0.25]);
set(gca,'FontSize',10);
subplot(1,2,2);imagesc(sftfalli);colorbar;ylim([0.5 3.5]);xlim([0.5 3.5])
yticks([1:1:3]);xticks([1:1:3]);
set(gca,'YTickLabel',{'1' '2','4'});set(gca,'Ydir','reverse');ylabel('TF');
set(gca,'XTickLabel',{'0.02' '0.08','0.16'});xlabel('SF');axis square;
hold on; title(['iviv' ' n=' num2str(sum(~isnan( od_out_iviv(:,9))))]);caxis([0 0.25]);
xtickangle(45);
set(gca,'FontSize',10);
 %% UMAP supplement
%% 
close all;
%UMAP projections Panel C and D
%load embedding
load('C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200523_embedding\reduced_data.mat')
plot_list = {'span','frac_vert'};
plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
%in vivo values non cyclic 
figure(1);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(1), 'Position', [200, 0, 200, 200])
figure(2);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(2), 'Position', [200, 0, 200, 200])
figure(3);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(3), 'Position', [200, 0, 200, 200])
figure(4);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(4), 'Position', [200, 0, 200, 200])
figure(5);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(5), 'Position', [200, 0, 200, 200])
figure(6);c=colorbar;caxis([0 16]);c.Ticks=[0:8:16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(6), 'Position', [200, 0, 200, 200])
figure(7);c=colorbar;caxis([0 1]);c.Ticks=[0:0.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(7), 'Position', [200, 0, 200, 200])
figure(9);c=colorbar;caxis([0 0.7]);c.Ticks=[0:0.35:0.7];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(9), 'Position', [200, 0, 200, 200])
figure(11);c=colorbar;caxis([0 0.6]);c.Ticks=[0:0.3:0.6];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(11), 'Position', [200, 0, 200, 200])
figure(12);c=colorbar;caxis([0 0.6]);c.Ticks=[0:0.3:0.6];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(12), 'Position', [200, 0, 200, 200])
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
%% Save
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Supp6\200524\span'
savepdf_SW(fn,1);
%% In vivo parameters
close all;
 plot_list = {'Sigmapref','Capeakpref','sad','pci','noise','SF','TF'};
 plotting_embedding_str(reduced_data, str, plot_list, 0,0, 'parula');
autoArrangeFigures;
%% 
figure(1);c=colorbar;caxis([2 4.3]);c.Ticks=[2:1:4.3];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(1), 'Position', [200, 0, 200, 200])
figure(2);c=colorbar;caxis([1 6]);c.Ticks=[1:2:6];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(2), 'Position', [200, 0, 200, 200])
 figure(3);c=colorbar;caxis([0 0.16]);c.Ticks=[0:0.08:0.16];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(3), 'Position', [200, 0, 200, 200])
 figure(4);c=colorbar;caxis([-4 1]);c.Ticks=[-4:2.5:1];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(4), 'Position', [200, 0, 200, 200])
figure(5);c=colorbar;caxis([-9.4 -2.4]);c.Ticks=[-9:3:-2.4];h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(5), 'Position', [200, 0, 200, 200])
figure(6);c=colorbar;caxis([0.02 0.16]);c.Ticks=[0.02:0.02:0.16]; h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(6), 'Position', [200, 0, 200, 200])
figure(7);c=colorbar;caxis([1 4]);c.Ticks=[1:1:4]; h = gca; h.XAxis.Visible = 'off';h = gca; h.YAxis.Visible = 'off';title('');set(gca,'FontSize',10);c.Label.FontSize=10;set(figure(7), 'Position', [200, 0, 200, 200])
%% Save
fn='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Figure4\200523_embedding\invivo'
savepdf_SW(fn,1); 
 
%% 
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
 %% 
 
 for i=1:length(str)
%whole map
tmp1=ex_map(3:end,:,i);
tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:));sum(tmp2(:))];
end
 %% Correlation matrix between input, morpho, pia

df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
com=[];com=[L23fr(morph_cells_id,1)  L23fr(morph_cells_id,2) L4fr(morph_cells_id,1)  L4fr(morph_cells_id,2)  span(morph_cells_id,1)...
    span(morph_cells_id,4) tot_input(morph_cells_id,1) tot_input(morph_cells_id,2) pia_input(morph_cells_id) morph_sel(morph_cells_id,:)]
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
     'Radial dis.','Total Length','Peak Nr. crossing','Dis. peak branch','Width'});set(gca,'FontSize',10)
 
 %% 
 
 com=[];com=[morph_sel(morph_cells_id,:) morph_parameters(morph_cells_id,[9]) morph_parameters(morph_cells_id,[10]) morph_parameters(morph_cells_id,[19]) morph_parameters(morph_cells_id,[20])]
 
 G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
 
 
;set(gca,'FontSize',10)
xticklabels({'AP Radial dis.','AP Total Length','AP Peak Nr. crossing','AP Dis. peak branch','AP Width/Height','AP Width','AP Height'...
     'BA Radial dis.','BA Total Length','BA Peak Nr. crossing','BA Dis. peak branch','BA Width/Height','BA Width','BA Height'});set(gca,'FontSize',10)
 xtickangle(45)
 yticklabels({'AP Radial dis.','AP Total Length','AP Peak Nr. crossing','AP Dis. peak branch','AP Width/Height','AP Width','AP Height'...
     'BA Radial dis.','BA Total Length','BA Peak Nr. crossing','BA Dis. peak branch','BA Width/Height','BA Width','BA Height'});set(gca,'FontSize',10)

 ytickangle(45)
 %% 
% corr_plot(com(:,2),com(:,3),[],{'','',''});xlabel('Apical Total Length','Color','k'); ylabel('Apical Peak number crossings','Color','k'); 
% corr_plot(com(:,2),com(:,6),[],{'','',''});xlabel('Apical Total Length','Color','k'); ylabel('Apical Width','Color','k'); 
 corr_plot(com(:,9),com(:,10),[],{'','',''});xlabel('Basal Total Length','Color','k'); ylabel('Basal Peak number crossings','Color','k'); 
 corr_plot(com(:,9),com(:,13),[],{'','',''});xlabel('Basal Total Length','Color','k'); ylabel('Basal Width','Color','k'); 
 %% 
 %% SF TF multiple barplots with Fraction 
par=morph_parameters(:,20)
g1=find(od_out_iviv(:,8)==0.02);
g2=find(od_out_iviv(:,8)==0.08);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
xticklabels({'0.02','0.08'});xtickangle(45);
%ylim([0 80]);yticks([0:40:80]);
ylabel('Extent Basal Y','Color','k'); ;set(gca,'FontSize',10); 
%% 

%% 
b=20
par1=[];
groups_idx=[];
g1=find(od_out_iviv(:,9)==1);
g2=find(od_out_iviv(:,9)==2);
g3=find(od_out_iviv(:,9)==4);
par1=vertcat(morph_parameters(g1,b),morph_parameters(g2,b),morph_parameters(g3,b))
groups_idx=vertcat(ones(length(g1),1)*1,ones(length(g2),1)*2,ones(length(g3),1)*3)
groups_idx(find(isnan(par1)))=[];
par1(find(isnan(par1)))=[];
[statsout] = barplot_sw(par1,groups_idx,{'','Extent Basal Y'});ylabel('L23 fraction')
%% 
par=span(:,3) 
g1=find(od_out_iviv(:,1)>0.25 &  od_out_iviv(:,4)>s2a & od_out_iviv(:,4)<s2b & od_out_iviv(:,8)==0.02);
g2=find(od_out_iviv(:,1)>0.25 &  od_out_iviv(:,4)>s2a & od_out_iviv(:,4)<s2b & od_out_iviv(:,8)==0.08);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
ylabel('Span','Color','k');xticklabels({'0.02','0.08'}) ;set(gca,'FontSize',10); 
%% 
par=L5fr(:,1) 
g1=find(od_out_iviv(:,1)>0.25 &  od_out_iviv(:,4)>s2a & od_out_iviv(:,4)<s2b & od_out_iviv(:,9)<=2);
g2=find(od_out_iviv(:,1)>0.25 &  od_out_iviv(:,4)>s2a & od_out_iviv(:,4)<s2b & od_out_iviv(:,9)==4);

[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
ylabel('fr','Color','k');xticklabels({'Slow','Fast'}) ;set(gca,'FontSize',10); xtickangle(45);
%% 

nvr=find([str(:).sftf_resp]==0 & [str(:).resp]==0);
nvra=find([str(:).resp]==0);
nvrb=find([str(:).sftf_resp]==1 & [str(:).resp]==0);
%nvr=find([str(:).resp]==0);
vr=find([str(:).resp]==1);
par=abs(ang_exL23(:,8)*69)
g1=vr;
g2=nvra;
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
ylabel('span','Color','k');xticklabels({'VR','NVR'}) ;set(gca,'FontSize',10); xtickangle(45);
%% Pia_all vs pia_input
od_id=[str(:).ODexpID]
pia_invivo=[str(:).pia_invivo]
figure;scatter(pia_invivo,od_id,20,pia_input,'filled');set(gcf,'color','w');ylabel('EXP ID');xlabel('Pia in vivo')
figure;histogram(pia_all);
%% 
iviv_id=find([str(:).iviv]==1)
u_id=unique(od_id(iviv_id))
p_i=pia_input(iviv_id)
p_ia=pia_invivo(iviv_id)
%% 

for i=1:31
temp=find(od_id(iviv_id)==u_id(i))
scale_di(i)=mean(p_ia(temp)-p_i(temp)')
end