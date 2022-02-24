%% L23 morpho/ephys structure
%filename=uipickfiles('FilterSpec','C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Stage4_allcells_structure')%pathname, you need uipickfiles function
filename=uipickfiles('FilterSpec','C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Stage5b_allinfo_structure')%pathname, you need uipickfiles function
load(char(filename));%load mat file

%% Read out which ones have morphology or not 
str_all=str;
for i=1:length(str_all)
    if isnan(str_all(i).RDA)==0
    morph(i)=1;    
    else
    morph(i)=0;
    end
end
morph_id=find(morph==1);
%% Read out morphology parameters for apical and basal tree
  data_morpho=[[str_all(morph_id).RDA]' [str_all(morph_id).LA]' [str_all(morph_id).PLA]' [[str_all(morph_id).BPA]'...
               [str_all(morph_id).BOA]' [str_all(morph_id).BLA]' [str_all(morph_id).WHA]' [str_all(morph_id).XSA]' [str_all(morph_id).YSA]' ...
               [str_all(morph_id).RDB]' [str_all(morph_id).LB]' [str_all(morph_id).PLB]' [str_all(morph_id).BPB]'...
               [str_all(morph_id).BOB]' [str_all(morph_id).BLB]' [str_all(morph_id).WHB]' [str_all(morph_id).XSB]' [str_all(morph_id).YSB]' [str_all(morph_id).NB]']];
%% Read out pial depth
for i=1:length(morph_id)
    temp=str_all(morph_id(i)).cellName
  name_morph{i}=temp(1:end-4);
if isnan(str_all(morph_id(i)).pialD) ==0    
pia_morpho(i)=str_all(morph_id(i)).pialD;
else
pia_morpho(i)=str_all(morph_id(i)).pia;   
end
end
length(unique(name_morph));
%% 
%%%%% Figure 1 %%%%%
%%  Histogram of pia distribution for morpho; Panel B
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_morpho,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
  ylabel('Pial depth (�m)');xlabel('Cell count');box off;
%  legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
% %     [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
 ylim([100 400])
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 50]);
% ;ylabel('Pial depth (�m)');xlabel('Cell count');box off;
%% sholl analyisis for example cell Panel C
[sex, ddex, sdex, XPex, YPex, ZPex, iDex] = sholl_tree(str_all(morph_id(86)).tr_basal, 20, '-s');
[sex1, ddex1, sdex1, XPex1, YPex1, ZPex1, iDex1] = sholl_tree(str_all(morph_id(86)).tr_apical, 20, '-s');hold on
%% Sholl Analysis APICAL for all cells Panel C
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 300, 225]);
for i=1:length(morph_id)
fig3=figure(3);   
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(str_all(morph_id(i)).tr_apical, 20, '-s');
max_a(i)=max(s);
temp=dd(find(s==max(s)));
dis_peaka(i)=temp(1);
close(fig3);
p1=plot(dd,s,'-k');
p1.Color(4) = 0.1;
hold on;
end
hold on;
% Sholl Analysis BASAL for all cells 
for i=1:length(morph_id)
fig3=figure(3); 
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(str_all(morph_id(i)).tr_basal, 20, '-s');
max_b(i)=max(s);
temp=dd(find(s==max(s)));
dis_peakb(i)=temp(1);
close(fig3);
p1=plot(dd,s,'-m');
p1.Color(4) = 0.1;
hold on;
end
box off;legend('Apical','Basal');legend boxoff;ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (�m)');set(gca,'FontSize',10);
%% highlight example cell Panel C
gca(fig1);hold on;p1=plot(ddex,sex,'-m','LineWidth',1);p1.Color(4) = 1;hold on;p1=plot(ddex1,sex1,'-k','LineWidth',1);p1.Color(4) = 1;hold on;
%% Correlation with pial depth: APICAL using sorted correlation plot, correct for multiple comparison, Panel D
stri={'Radial Dis._{max} (�m)','Total Length (�m)','Path Length_{max} (�m)','Branch Points','Branch Order_{max}','Branch Length_{mean} (�m)','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch (�m)'};
sorted_correlation([data_morpho(:,1:9) max_a' dis_peaka'],pia_morpho',stri,'k',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Apical');
%% Correlation with pial depth: BASAL using sorted correlation plot, panel E
stri={'Radial Dis._{max} (�m)','Total Length (�m)','Path Length_{max} (�m)','Branch Points','Branch Order_{max}','Branch Length_{mean} (�m)','Width / Height','Width','Height','Nr. Branches','Peak Nr. cross','Dis. peak branch (�m)'}
sorted_correlation([data_morpho(:,10:19) max_b' dis_peakb'],pia_morpho',stri,'m',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Basal');
%% PCA for apical morphology
[coeff_morph_a,score_morph_a,latent_morph_a,~,explained_morph_a,mu] = pca(zscore([data_morpho(:,1:9) max_a' dis_peaka']));
var_exp(explained_morph_a,[],[]);
%show only first three as suggested from TR
hold on;xlim([0.5 3.5]);
%% PCA for basal morphology
[coeff_morph_b,score_morph_b,latent_morph_b,~,explained_morph_b,mu] = pca(zscore([data_morpho(:,10:19) max_b' dis_peakb']));
var_exp(explained_morph_b,[],[]); 

%% PCA for both apical and basal
[coeff_morph,score_morph,latent_morph,~,explained_morph,mua] = pca(zscore([data_morpho max_a' dis_peaka' max_b' dis_peakb']));
var_exp(explained_morph_a,[],[]);
%% Plot the 1 corrleation each for Apical raw and PCA Panel F and G
tr = rescale(pia_morpho);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 500]);
subplot(2,1,1);scatter(data_morpho(:,7)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');title('Apical');
ylabel('Relative pial position');xlabel('Width / Height');
text(4,0.99,['r=-0.64']);text(4,0.9,['p<0.001']);set(gca,'FontSize',10);xlim([0 4.5]);ylim([0 1])
 P = polyfit(data_morpho(:,7)',tr,1);
    yfit = P(1)*data_morpho(:,7)'+P(2);
    hold on;
    plot(data_morpho(:,7)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(2,1,2);
scatter(score_morph_a(:,2),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
xlim([-4 5]);
ylabel('Relative pial position');xlabel('PC2 (25.95%)');
%xlim([25 85]);xticks([25:20:85]);ref.XData=[25 85];
text(-3.5,0.99,['r=0.8']);text(-3.5,0.9,['p<0.001']);set(gca,'FontSize',10);P=[];yfit=[];
P = polyfit(score_morph_a(:,2)',tr,1);
yfit = P(1)*score_morph_a(:,2)'+P(2);hold on;
plot(score_morph_a(:,2),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
 %% Numbers for correlations 
 %apical
 [r p]=corrcoef(score_morph_a(:,4),tr);
 ra=r(2)
 pa=p(2)
 %basal
 [r p]=corrcoef(score_morph_b(:,3),tr);
 rb=r(2)
 pb=p(2)
 %% Ephys Figure 2
%load example traces for Panel A
load('D:\Munich_backup_harddriveNEW\THESIS\Ephys\ephys.mat')
%example cells
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 800, 175]);
subplot(1,3,1)
%Rheobase +30 pA
%up: SW1504080002:151.8
plot(ephys.spike_traces{:,:,48}(:,9),'k');box off;ylim([-80 80]);;yticks([-80:40:80]);
%mid: SW150403004: 222.2
subplot(1,3,2)
plot(ephys.spike_traces{:,:,46}(:,14),'k');box off;ylim([-80 80]);;yticks([-80:40:80]);
%low: SW150403004: 342
subplot(1,3,3)
 plot(ephys.spike_traces{:,:,34}(:,17),'k');ylim([-80 80]);yticks([-80:40:80]);
 %ephys.spike_traces{:,:,34}(1,17); 
 box off;
%% %% Read out which ones have ephys or not 
ephys=[];
str_all=str;
for i=1:length(str_all)
    if isnan(str_all(i).Vmin)==0
    ephys(i)=1;    
    else
    ephys(i)=0;
    end
end
ephys_id=find(ephys==1);
for i=1:length(ephys_id)
   temp=str_all(ephys_id(i)).cellName
  name_ephys{i}=temp(1:end-4);
end
length(unique(name_ephys));
%% Read out ephys parameters for passive and active
  data_ephys=[[str_all(ephys_id).Vmin]' [str_all(ephys_id).Vpeak]' [str_all(ephys_id).Vthresh]' [[str_all(ephys_id).Vslope]'...
               [str_all(ephys_id).Vhalf]' [str_all(ephys_id).Vamp]' [str_all(ephys_id).AHP]' [str_all(ephys_id).APrise]' [str_all(ephys_id).APfall]' ...
               [str_all(ephys_id).APbwidth]' [str_all(ephys_id).APhwidth]' [str_all(ephys_id).APlate]' [str_all(ephys_id).Vrest]'...
               [str_all(ephys_id).tau]' [str_all(ephys_id).Rin]' [str_all(ephys_id).Sag]' [str_all(ephys_id).Rheo]' [str_all(ephys_id).APfreq]']];
%% %% Read out pial depth for ephys
for i=1:length(ephys_id)
pia_ephys(i)=str_all(ephys_id(i)).pia;    
end
%%  Histogram of pia distribution for ephys; Panel B
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_ephys,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
 ylabel('Pial depth (�m)');xlabel('Cell count');box off;ylim([100 400]);
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 40]);
 %% Correlation with pial depth: Passive using sorted correlation plot;multiple comparison corrected; Panel C
stri={'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'}
sorted_correlation([data_ephys(:,13:17)],pia_ephys',stri,'k',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Passive');
%% Correlation with pial depth: Active using sorted correlation plot;multiple comparison corrected; Panel D
stri={'APV_{min} (mV)','APV_{peak} (mV)','APV_{thresh} (mV)','APV_{slope} (mV)','APV_{half} (mV)','APV_{amp} (mV)','AHP (mV)','APfreq (Hz)'}
sorted_correlation([data_ephys(:,[1:7 18])],pia_ephys',stri,[0 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Active');
%% PCA for passive ephys
ephys_sel=[];ephys_sel=[data_ephys(:,13:17)];score_ephys=[];
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
var_exp(explained_ephys,[],[]); 
%% 
ephys_sel=[];ephys_sel=[data_ephys(:,[1:7 18])];score_ephys=[];
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
var_exp(explained_ephys,[],[]);

%% Plot the 1 corrleation each for passive and active
tr=[];tr = rescale(pia_ephys);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_ephys(:,14)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
title('Passive');text(55,0.99,'r=-0.29');text(55,0.9,'p<0.001');hold on;
xlim([10 80]);xticks([0:20:80]);ylabel('Relative pial position');xlabel('Tau (ms)')
P=[];yfit=[];
 P = polyfit(data_ephys(:,14)',tr,1);
    yfit = P(1)*data_ephys(:,14)'+P(2);
    hold on;
    plot(data_ephys(:,14)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(1,2,2);
scatter(score_ephys(:,1),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');hold on;xticks([-5:2:5]);
ylabel('Relative pial position');xlabel('PC1 (53%)')
xlim([-5 5]);%xticks([40:80:280]);
text(3.5,0.99,'r=-0.32');text(3.5,0.9,'p<0.001');
P=[];yfit=[];
 P = polyfit(score_ephys(:,1)',tr,1);
    yfit = P(1)*score_ephys(:,1)'+P(2);
    hold on;
    plot(score_ephys(:,1)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
 %% MORPHO AND EPHYS
%% Overlapping morpho ephys
[~,ia,ib] =intersect(ephys_id,morph_id)
%% Read out 33 cells that overlap
data_ephys_both=[];
data_morph_both=[];
data_morph_both=[data_morpho(ib,:)];
data_ephys_both=[data_ephys(ia,:)];
pia_both=[pia_morpho(ib)];
trm=morph_id(ib);
for i=1:length(trm)
    temp=str_all(trm(i)).cellName
  name_mephys{i}=temp(1:end-4);
end
length(unique(name_mephys))   
%% Correlation of both
%Passive with apical parameters
ap_morph=[data_morph_both(:,1:9) max_a(ib)' dis_peaka(ib)'];
G2=correlation_matrix([ap_morph data_ephys_both(:,[13:17])],0);
%% Apical
ap_morph=[];
ap_morph=[data_morph_both(:,1:9) max_a(ib)' dis_peaka(ib)'];
%ap_morph=[data_morph_both(:,10:end) max_b(ib)' dis_peakb(ib)'];
G2=correlation_matrix([ap_morph data_ephys_both(:,[1:7 18])],0);
% %% Basal
% ap_morph=[];
% ap_morph=[data_morph_both(:,10:end) max_b(ib)' dis_peakb(ib)'];
% G2=correlation_matrix([ap_morph data_ephys_both(:,[13:17])],0);
%% correlation matrix; non multiple comparison corrected; Panel G
[R2,P2]=corrcoef([ap_morph data_ephys_both(:,[13:17])],'rows','pairwise');
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 300])
imagesc(R2(12:end,1:11));c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');colormap(cmap);caxis([-1 1]);xticks([1:1:11]);yticks([1:1:5])
yticklabels({'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'});set(gca,'FontSize',10);
xticklabels({'Radial Dis._{max}','Total Length','Path Length_{max}','Branch Points','Branch Order_{max}','Branch Length_{mean}','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch'});
xtickangle(45);set(gca,'FontSize',10);
%% Input connectivity Figure 3 
%% Input maps
str_all=str;
for i=1:length(str)
    if isnan(str(i).Cluster_id)==0
    imaps(i)=1;    
    else
    imaps(i)=0;
    end
end
imaps_id=find(imaps==1);
imaps_id(44)=[];
for i=1:length(imaps_id)
    temp=str_all(imaps_id(i)).cellName
  name_maps{i}=temp(1:end-4);
end
length(unique(name_maps));  
%% read out 16x16 maps 
for i=1:length(imaps_id)
ex_map(:,:,i) = str(imaps_id(i)).subpixel_excMap;
in_map(:,:,i) = str(imaps_id(i)).subpixel_inhMap;
frh(i,:)=str(imaps_id(i)).frac_horz;
pia_input(i,:)=[str(imaps_id(i)).pialD]';
cellid_inp(i)=[str(imaps_id(i)).cellID]';
end
for i=1:length(imaps_id)
ex_map_raw(:,:,i) = str(imaps_id(i)).subpixel_raw_excMap;
in_map_raw(:,:,i) = str(imaps_id(i)).subpixel_raw_inhMap;
end
%Calculate simple difference between maps
diff_map=ex_map-in_map;
str_imap=str(imaps_id);
%% calculate fraction vertical and horizontal
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv L23h L4h L5h pialD layer_assign] = iviv_profiles(1:147,str_imap); combine
for i=1:length(imaps_id)
frv(i,:)=[[0 0] frac_exv(i,:) frac_inv(i,:)];
end
%% calculate fraction per layer + span
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
[span(:,1),span(:,2),span(:,3),span(:,4),span(:,5),span(:,6)] = span_perLayer(ex_map,in_map,1:147);
spandL23=span(:,1)-span(:,4);
spandL4=span(:,2)-span(:,5);
spandL5=span(:,3)-span(:,6);
%% Show three map examples, Panel A
% cell 108,  165 microm depth
plot_avg_maps(str_imap,108,ex_map_raw,in_map_raw,pia_input,1,0,[]);
% cell 110, 180 microm depth
plot_avg_maps(str_imap,110,ex_map_raw,in_map_raw,pia_input,1,0,[]);
% cell 80, 320 microm depth
plot_avg_maps(str_imap,80,ex_map_raw,in_map_raw,pia_input,1,0,[]);
%%  Histogram of pia distribution for inputs, Panel B
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_input,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
 ylabel('Pial depth (�m)');xlabel('Cell count');box off;ylim([100 400])
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 40]);  
 %% Correlation with pial depth: Vertical fraction; Panel C, multiple comparison corrected
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
%sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2),diffL23fr,diffL4fr,diffL4fr],pia_input,stri,[0.5 0.5 0.5],1)
sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2)],pia_input,stri,[0.5 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Vertical fraction');
%% Correlation with pial depth: Maximum horizontal extent; Panel D, multiple comparison corrected
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
 %sorted_correlation([span,spandL23,spandL4,spandL5],pia_input,stri,[0.5 0.5 0.5],0)
sorted_correlation([span],pia_input,stri,[0.5 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Maximal horizontal extent');
%% Plot most significant corrleations vertical fractions; Panel E
tr=[];tr = rescale(pia_input);
%L4 EX
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
subplot(3,1,1);
scatter(L4fr(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r'); xlim([0 0.7])
text(0,0.99,'r=0.41');text(0,0.9,'p<0.001');hold on;
xticks([0:0.15:0.7]);xlim([-0.05 0.7])
P=[];yfit=[];
P = polyfit(L4fr(:,1),tr,1);
    yfit = P(1)*L4fr(:,1)+P(2);
    hold on;
    plot(L4fr(:,1),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
    
%L4 IN
subplot(3,1,2);
scatter(L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');xlim([-0.05 0.6]);
text(0,0.99,'r=0.3');text(0,0.9,'p<0.001');xlabel('Fraction');xticks([0:0.15:0.6]);
P=[];yfit=[];
P = polyfit(L4fr(:,2),tr,1);
    yfit = P(1)*L4fr(:,2)+P(2);
    hold on;
    plot(L4fr(:,2),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); ylabel('Relative pial position');

%L4 EX-IN
subplot(3,1,3);
scatter(L4fr(:,1)'-L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 xlim([-0.2 0.6]);
 xlabel('\Delta EX-IN vertical fraction');
 hold on;line([0 0], [0 1],'Color','k','LineStyle','--');xlim([-0.4 0.5]);
%  hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
% ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Fraction')
% ref.XData=[0 0.6];xticks([0:0.15:0.6]);
 text(-0.3,0.99,'r=0.15');text(-0.3,0.9,'n.s.');
 %% Plot most significant corrleations horizontal extent; Panel F
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
%L2/3 EX
subplot(3,1,1);
scatter(span(:,1)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');xlim([0 15*69]);text(800,0.99,'r=-0.3');text(800,0.9,'p<0.001');
xticks([0:500:1000]);
P=[];yfit=[];
P = polyfit(span(:,1)*69,tr,1);
    yfit = P(1)*span(:,1)*69+P(2);
    hold on;
    plot(span(:,1)*69,yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
%L2/3 IN
subplot(3,1,2);
scatter(span(:,4)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');set(gca,'FontSize',10);xlim([0 15*69]);
text(800,0.99,'r=-0.22');text(800,0.9,'n.s.');yticks([0:0.25:1]);ylabel('Relative pial position');xlabel('Horizontal Extent (�m)')
xticks([0:500:1000]);

%L2/3 EX-IN
subplot(3,1,3);
scatter(spandL23*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
xlim([-500 500]);
xlabel('\Delta EX-IN vertical fraction');

xlabel('\Delta EX-IN horizontal extent')
hold on;line([0 0], [0 1],'Color','k','LineStyle','--');set(gca,'FontSize',10)
text(400,0.99,'r=-0.16');text(400,0.9,'n.s.');
%% PCA with whole maps; Panel G
[coeff_ex,score_ex,coeff_in,score_in,coeff_com,score_com] = map_align_PCA(str,imaps_id); 
%% Plot correlation matrix of map scores; Panel H
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
scatter(score_com(:,1),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('Relative pial position');xlabel('PC1_{com}');set(gca,'FontSize',10);
yticks([0:0.25:1]);
text(-2,0.99,'r=0.44');text(-2,0.9,'p<0.001');
P=[];yfit=[];
P = polyfit(score_com(:,1),tr,1);
    yfit = P(1)*score_com(:,1)+P(2);
    hold on;
    plot(score_com(:,1),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
%% Input maps to morphology
%% Overlapping morpho input
ia=[];ib=[];
[~,ia,ib] =intersect(imaps_id,morph_id)
%% Read out cells that overlap morpho input
data_input=[L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2),diffL23fr,diffL4fr,diffL4fr...
    span,spandL23,spandL4,spandL5];
data_mi_morph=[];data_mi_input=[];
data_mi_morph=[data_morpho(ib,:)];
data_mi_input=[data_input(ia,:)];
pia_mi=[pia_morpho(ib)];
trmm=morph_id(ib);
for i=1:length(trmm)
    temp=str_all(trmm(i)).cellName
  name_mmaps{i}=temp(1:end-4);
end
length(unique(name_mmaps))  
%% matrix
com=[];com=[data_mi_morph data_mi_input]
G2=correlation_matrix([data_mi_morph data_mi_input],0);

%% Basal
ba_input=[data_mi_morph(:,10:19) max_b(ib)' dis_peakb(ib)'];
[R2,P2]=corrcoef([ba_input data_mi_input],'rows','pairwise');
G2=correlation_matrix([ba_input data_mi_input],0);
%% Apical
ap_input=[data_mi_morph(:,1:9) max_a(ib)' dis_peaka(ib)'];
[R2,P2]=corrcoef([ap_input data_mi_input],'rows','pairwise');
G2=correlation_matrix([ap_input data_mi_input],0);
%% combined
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 300])
imagesc(R2(13:end,1:12));c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:11]);yticks([1:1:5])

%% Plot example correlation for Panel I
tr=[]; tr=rescale(pia_mi);
cmap='plasma'
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 300])
scatter(data_mi_input(:,10)*69,ba_input(:,7),30,tr,'filled');
ylabel('Basal Width / Height');xlabel('Input Horizontal extent L2/3_{EX} (�m)');set(gca,'FontSize',10);
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[0 2.6];;ref.Color='k'; 
 ylim([0.5 2.6]);xlim([0 1100]);
 colormap(cmap);box off; hold on; c=colorbar;%c.Label.String = 'Relative pial position'
 text(100,2.5,'r=0.34');text(100,2.4,'p<0.001');
   P=[];yfit=[];
 P = polyfit(data_mi_input(:,10)'*69,ba_input(:,7)',1);
    yfit = P(1)*data_mi_input(:,10)'*69+P(2);
    hold on;
    plot(data_mi_input(:,10)'*69,yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w');   %% 
 
%% INVIVO, Figure 4

%% in vivo in vitro comparison
% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
directory=out_dir;% use cobined date structure named 
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
[od_out sftf_out sftf_sel sftf_pref spon_out pia_all delta_Ca fit_Ca  Ca_peak_od] = concat_invivo(L23_PC_new);         
%% calculate R2 of fit
oris=[0:45:315];
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
 %% 
%% Sorted correlation checked with multiple comparison; Panel C left
tr=[];tr=rescale(pia_all);
stri={'gOSI','gDSI','ODI','Tuning Width','R / R0_{max}'};
a=[];a=find(od_out(:,1)>0.25);
a=[];a=r2_all'<0.3;
gOSI=od_out(:,1);
gOSI(a)=NaN;
a=[];a=r2_all'<0.3;
gDSI=od_out(:,2);
gDSI(a)=NaN;
a=[];a=r2_all'<0.3;
TW=od_out(:,7);
TW(a)=NaN;
sorted_correlation([gOSI gDSI od_out(:,3) TW max(delta_Ca(:,:),[],2)],tr',stri,[0.5 0.5 0.5],1);
xlabel('Correlation with pial depth');
%% perform PCA remove NaNs;
invivo_feat=[gOSI gDSI od_out(:,3) TW max(delta_Ca(:,:),[],2) tr']
invivo_feat(any(isnan(invivo_feat),2),:) = [];
[coeff_invivo,score_invivo,latent_invivo,~,explained_invivo,mu] = pca(zscore(invivo_feat(:,1:5)));
%% Panel C right
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 300])
imagesc(coeff_invivo([3 5 4 1 2],1:3));
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:3]);yticks([1:1:5]);colorbar;
xticklabels({'PC1','PC2','PC3'});
yticklabels({'ODI','R / R0_{max}','Tuning Width','gOSI','gDSI'});
set(gca,'FontSize',10);
%% Show PC2 vs pial depth; Panel D
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 300])
scatter(score_invivo(:,2),invivo_feat(:,6),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');hold on;
ylabel('Relative pial position');xlabel('PC2 (22%)')
text(4,0.99,'r=0.1');text(4,0.9,'p<0.01');axis square
P=[];yfit=[];
 P = polyfit(score_invivo(:,2),invivo_feat(:,6),1);
    yfit = P(1)*score_invivo(:,2)'+P(2);
    hold on;
    plot(score_invivo(:,2)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
    %% Responiveness per depth in vivo all using sum, split L2/3 in half
g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.5);
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.5);  
par=[];
par=max(delta_Ca(:,:),[],2);
[statsout]=dual_barplot(par,g1,g3,2);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Violin boxplot;  split L2/3 in half: Responsiveness amplitude
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 260, 400]);
hold on;violin({par(g1) par(g3)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','m'); box off; ylabel('dR / R0_{max}'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g3));
xlabel('Relative pial position');hold on;
%% 