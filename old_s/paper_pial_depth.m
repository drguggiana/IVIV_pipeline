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
           %% 
%            
%        mean_ap= nanmean([data_morpho(:,1:9) max_a' dis_peaka']);
%          mean_ba= nanmean([data_morpho(:,10:19) max_b' dis_peakb']);
%         
%           
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
length(unique(name_morph))
%%  Histogram of pia distribution for morpho
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_morpho,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
  ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%  legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
% %     [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
 ylim([100 400])
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 50]);
% ;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%% Sholl Analysis APICAL
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
% Sholl Analysis BASAL
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
box off;legend('Apical','Basal');legend boxoff;ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%% 

[sex, ddex, sdex, XPex, YPex, ZPex, iDex] = sholl_tree(str_all(morph_id(86)).tr_basal, 20, '-s');
[sex1, ddex1, sdex1, XPex1, YPex1, ZPex1, iDex1] = sholl_tree(str_all(morph_id(86)).tr_apical, 20, '-s');
%% 

hold on
p1=plot(ddex,sex,'-m','LineWidth',1);
p1.Color(4) = 1;
hold on;
p1=plot(ddex1,sex1,'-k','LineWidth',1);
p1.Color(4) = 1;
hold on;
%% Correlation with pial depth: APICAL using sorted correlation plot
stri={'Radial Dis._{max} (µm)','Total Length (µm)','Path Length_{max} (µm)','Branch Points','Branch Order_{max}','Branch Length_{mean} (µm)','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch (µm)'};
sorted_correlation([data_morpho(:,1:9) max_a' dis_peaka'],pia_morpho',stri,'k',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Apical');
%% Correlation with pial depth: BASAL using sorted correlation plot
stri={'Radial Dis._{max} (µm)','Total Length (µm)','Path Length_{max} (µm)','Branch Points','Branch Order_{max}','Branch Length_{mean} (µm)','Width / Height','Width','Height','Nr. Branches','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,10:19) max_b' dis_peakb'],pia_morpho',stri,'m',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Basal')
% %% Sliding winddow analysis APICAL
% [p_wha1 p_wha2] = sliding_window_comp(pia_morpho,data_morpho(:,3));
% %[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,max_a);
% %% Sliding winddow analysis BASAL
% [p_blb1 p_blb2] = sliding_window_comp(pia_morpho,data_morpho(:,15));
% [p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,data_morpho(:,13));
%% PCA for apical morphology
[coeff_morph_a,score_morph_a,latent_morph_a,~,explained_morph_a,mu] = pca(zscore([data_morpho(:,1:9) max_a' dis_peaka']));
var_exp(explained_morph_a,[],[]); 
%% 
%% PCA for basal morphology
[coeff_morph_b,score_morph_b,latent_morph_b,~,explained_morph_b,mu] = pca(zscore([data_morpho(:,10:19) max_b' dis_peakb']));
var_exp(explained_morph_b,[],[]); 

%% Plot the 1 corrleation each for Apical and BAsal
tr = rescale(pia_morpho);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 500]);
subplot(2,1,1);
scatter(data_morpho(:,7)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 title('Apical');
%  hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[0 4];xticks([0:1:4]);
% ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';
ylabel('Relative pial position');xlabel('Width / Height');
text(4,0.99,['r=-0.64']);text(4,0.9,['p<0.001']);set(gca,'FontSize',10);xlim([0 4.5]);ylim([0 1])
 P = polyfit(data_morpho(:,7)',tr,1);
    yfit = P(1)*data_morpho(:,7)'+P(2);
    hold on;
    plot(data_morpho(:,7)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(2,1,2);
scatter(score_morph_a(:,2),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
% hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ;ref.XData=[-4 5];%xticks([0.25:1:4]);
xlim([-4 5]);
%ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';
ylabel('Relative pial position');xlabel('PC2 (25.95%)');
%xlim([25 85]);xticks([25:20:85]);ref.XData=[25 85];
text(-3.5,0.99,['r=0.8']);text(-3.5,0.9,['p<0.001']);set(gca,'FontSize',10);
P=[];yfit=[];
 P = polyfit(score_morph_a(:,2)',tr,1);
    yfit = P(1)*score_morph_a(:,2)'+P(2);
    hold on;
    plot(score_morph_a(:,2),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

%% Ephys
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
 %ephys.spike_traces{:,:,34}(1,17)
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
length(unique(name_ephys))
%% Read out morphology parameters for apical and basal tree
  data_ephys=[[str_all(ephys_id).Vmin]' [str_all(ephys_id).Vpeak]' [str_all(ephys_id).Vthresh]' [[str_all(ephys_id).Vslope]'...
               [str_all(ephys_id).Vhalf]' [str_all(ephys_id).Vamp]' [str_all(ephys_id).AHP]' [str_all(ephys_id).APrise]' [str_all(ephys_id).APfall]' ...
               [str_all(ephys_id).APbwidth]' [str_all(ephys_id).APhwidth]' [str_all(ephys_id).APlate]' [str_all(ephys_id).Vrest]'...
               [str_all(ephys_id).tau]' [str_all(ephys_id).Rin]' [str_all(ephys_id).Sag]' [str_all(ephys_id).Rheo]' [str_all(ephys_id).APfreq]']];

%% %% Read out pial depth
for i=1:length(ephys_id)
pia_ephys(i)=str_all(ephys_id(i)).pia;    
end
%%  Histogram of pia distribution for ephys
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_ephys,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
  ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%  legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
% %     [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
 ylim([100 400])
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 40]);
% ;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%% Correlation with pial depth: Passive using sorted correlation plot
stri={'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'}
sorted_correlation([data_ephys(:,13:17)],pia_ephys',stri,'k',1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Passive');
%% Correlation with pial depth: Active using sorted correlation plot
stri={'APV_{min} (mV)','APV_{peak} (mV)','APV_{thresh} (mV)','APV_{slope} (mV)','APV_{half} (mV)','APV_{amp} (mV)','AHP (mV)','APfreq (Hz)'}
sorted_correlation([data_ephys(:,[1:7 18])],pia_ephys',stri,[0 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Active');
% %% Sliding winddow analysis Passive
% [p_wha1 p_wha2] = sliding_window_comp(pia_ephys,data_ephys(:,17));
% %% Sliding winddow analysis Active
% [p_blb1 p_blb2] = sliding_window_comp(pia_ephys,data_ephys(:,4));
%% PCA for passive ephys
ephys_sel=[];ephys_sel=[data_ephys(:,13:17)];
score_ephys=[];
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
var_exp(explained_ephys,[],[]); 
%% Plot the 1 corrleation each for passive and active
tr=[];tr = rescale(pia_ephys);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_ephys(:,14)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 title('Passive');text(55,0.99,'r=-0.29');text(55,0.9,'p<0.001');
 hold on;
 %ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];
 xlim([10 80]);xticks([0:20:80]);
%ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';
ylabel('Relative pial position');xlabel('Tau (ms)')

%ref.XData=[15 70];xticks([10:20:80]);
P=[];yfit=[];
 P = polyfit(data_ephys(:,14)',tr,1);
    yfit = P(1)*data_ephys(:,14)'+P(2);
    hold on;
    plot(data_ephys(:,14)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(1,2,2);
scatter(score_ephys(:,1),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 hold on;
 %ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[30 90];
 xticks([-5:2:5]);
%ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';
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
%% Basal

ap_morph=[];
ap_morph=[data_morph_both(:,10:end) max_b(ib)' dis_peakb(ib)'];
G2=correlation_matrix([ap_morph data_ephys_both(:,[13:17])],0);
%% Fugyre correlation matrix

[R2,P2]=corrcoef([ap_morph data_ephys_both(:,[13:17])],'rows','pairwise');

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 300])
imagesc(R2(12:end,1:11));c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:11]);yticks([1:1:5])
yticklabels({'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'});set(gca,'FontSize',10);
xticklabels({'Radial Dis._{max}','Total Length','Path Length_{max}','Branch Points','Branch Order_{max}','Branch Length_{mean}','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch'});
xtickangle(45);set(gca,'FontSize',10);


%% PCA for morpho ephys
% 
% ephys_sel=data_ephys(:,[13 14 15 17]);
% morph_sel=[data_morpho(:,[1 3 5 6 7 8 9]) max_a' dis_peaka'];
% [coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
% [coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(zscore(morph_sel));
% var_exp(explained_ephys,[],[]); 
% var_exp(explained_morph,[],[]); 
%% 
 tr=[];
 tr = rescale(pia_ephys);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 200]);
scatter(score_ephys(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor',[0 0.5 0.5]);
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 xlim([-6 5]);ref.XData=[-6 5];xticks([-6:2:5]);ref.YData=[1 0];
 ref.Color='k';ylabel('Relative pial position');xlabel('PC1_{e}');
 tr=[];
 tr = rescale(pia_morpho);
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 200]);
scatter(score_morph(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor',[0.5 0.5 0.5]);
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 xlim([-6 5]);ref.XData=[-6 5];xticks([-6:2:5]);ref.YData=[0 1];
 ref.Color='k';ylabel('Relative pial position');xlabel('PC1_{e}');
%% Input maps

%% 

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
length(unique(name_maps))  
%% 
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
%% 
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv L23h L4h L5h pialD layer_assign] = iviv_profiles(1:147,str_imap);
%% 
for i=1:length(imaps_id)
frv(i,:)=[[0 0] frac_exv(i,:) frac_inv(i,:)];
end
%% 

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
spandL23=span(:,1)-span(:,4);
spandL4=span(:,2)-span(:,5);
spandL5=span(:,3)-span(:,6);
%%  Histogram of pia distribution for inputs
 close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 225]);
 hold on;histogram(pia_input,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
  ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%  legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
% %     [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
 ylim([100 400])
 set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 40]);
%% Show three map examples
% cell 180
plot_avg_maps(str_imap,110,ex_map_raw,in_map_raw,pia_input,1,0,[]);
%% 165
plot_avg_maps(str_imap,108,ex_map_raw,in_map_raw,pia_input,1,0,[]);
%% 320
plot_avg_maps(str_imap,80,ex_map_raw,in_map_raw,pia_input,1,0,[]);
%% %Pial depth correlation with EX IN L23, L4, L5 for panel F
plot_fraction(L23fr,L4fr,L5fr,pia_input);
%% 
%% Correlation with pial depth: Vertical fraction
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
%sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2),diffL23fr,diffL4fr,diffL4fr],pia_input,stri,[0.5 0.5 0.5],1)
sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2)],pia_input,stri,[0.5 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Vertical fraction');
%% Correlation with pial depth: Maximum horizontal extent
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
 %sorted_correlation([span,spandL23,spandL4,spandL5],pia_input,stri,[0.5 0.5 0.5],0)
sorted_correlation([span],pia_input,stri,[0.5 0.5 0.5],1)
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Maximal horizontal extent');
%% Plot most significant corrleations vertical fractions
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
subplot(3,1,1);
scatter(L4fr(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');
 %set(gca, 'YDir','reverse');
 xlim([0 0.7])
 text(0,0.99,'r=0.41');text(0,0.9,'p<0.001');
 hold on;
 %ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];
%ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';

%ref.XData=[0 0.7];
xticks([0:0.15:0.7]);xlim([-0.05 0.7])

P=[];yfit=[];
 P = polyfit(L4fr(:,1),tr,1);
    yfit = P(1)*L4fr(:,1)+P(2);
    hold on;
    plot(L4fr(:,1),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(3,1,2);
scatter(L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
 %set(gca, 'YDir','reverse');
xlim([-0.05 0.6]);
 text(0,0.99,'r=0.3');text(0,0.9,'p<0.001');
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
%ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';
xlabel('Fraction');
%ref.XData=[0 0.6];
xticks([0:0.15:0.6]);

P=[];yfit=[];
 P = polyfit(L4fr(:,2),tr,1);
    yfit = P(1)*L4fr(:,2)+P(2);
    hold on;
    plot(L4fr(:,2),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
ylabel('Relative pial position');

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
%% Plot most significant corrleations horizontal extent
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
subplot(3,1,1);
scatter(span(:,1)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');
 %set(gca, 'YDir','reverse');
 xlim([0 15*69]);
 text(800,0.99,'r=-0.3');text(800,0.9,'p<0.001');
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';
 %ylabel('Relative pial position');xlabel('Horizontal Extent (µm)');set(gca,'FontSize',10)
 %ref.XData=[0 15*69];
 xticks([0:500:1000]);
 
 P=[];yfit=[];
 P = polyfit(span(:,1)*69,tr,1);
    yfit = P(1)*span(:,1)*69+P(2);
    hold on;
    plot(span(:,1)*69,yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

subplot(3,1,2);
scatter(span(:,4)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');set(gca,'FontSize',10)
 %set(gca, 'YDir','reverse');
 xlim([0 15*69]);
 text(800,0.99,'r=-0.22');text(800,0.9,'n.s.');
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
% ref.YData=[1 0];
yticks([0:0.25:1]);
ylabel('Relative pial position');xlabel('Horizontal Extent (µm)')
 %ref.XData=[0 15*69];
 xticks([0:500:1000]);

subplot(3,1,3);
scatter(spandL23*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 xlim([-500 500]);
 xlabel('\Delta EX-IN vertical fraction');

 %ref.XData=[0 4];xticks([0:1:4]);
%ylabel('Relative pial position');
xlabel('\Delta EX-IN horizontal extent')
 hold on;line([0 0], [0 1],'Color','k','LineStyle','--');set(gca,'FontSize',10)
  text(400,0.99,'r=-0.16');text(400,0.9,'n.s.');
%% PCA with whole maps
 [coeff_ex,score_ex,coeff_in,score_in,coeff_com,score_com] = map_align_PCA(str,imaps_id);
%% Plot correlation matrix of map scores
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
scatter(score_com(:,1),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('Relative pial position');xlabel('PC1_{com}');set(gca,'FontSize',10);
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[0 1];
 yticks([0:0.25:1]);
 %ref.Color='k'; 
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
%% 
com=[];com=[data_mi_morph data_mi_input]
G2=correlation_matrix([data_mi_morph data_mi_input],0);

%% Basal
ba_input=[data_mi_morph(:,10:19) max_b(ib)' dis_peakb(ib)'];
[R2,P2]=corrcoef([ba_input data_mi_input],'rows','pairwise');
G2=correlation_matrix([ba_input data_mi_input],0);
%% 
ap_input=[data_mi_morph(:,1:9) max_a(ib)' dis_peaka(ib)'];
[R2,P2]=corrcoef([ap_input data_mi_input],'rows','pairwise');
G2=correlation_matrix([ap_input data_mi_input],0);
%% 

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 300])
imagesc(R2(13:end,1:12));c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:11]);yticks([1:1:5])

%% Plot specific input and morpho
tr=[]; tr=rescale(pia_mi);
corr_plot(data_mi_input(:,10),ba_input(:,8),tr,{'Input Horizontal extent L2/3_{EX}','Width / Height_{basal}','Relative pial position'});
%% 
tr=[]; tr=rescale(pia_mi);
cmap='plasma'
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 300])
scatter(data_mi_input(:,10)*69,ba_input(:,7),30,tr,'filled');
ylabel('Basal Width / Height');xlabel('Input Horizontal extent L2/3_{EX} (µm)');set(gca,'FontSize',10);
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[0 2.6];;ref.Color='k'; 
 ylim([0.5 2.6]);xlim([0 1100]);
 colormap(cmap);box off; hold on; c=colorbar;%c.Label.String = 'Relative pial position'

 text(100,2.5,'r=0.34');text(100,2.4,'p<0.001');
   P=[];yfit=[];
 P = polyfit(data_mi_input(:,10)'*69,ba_input(:,7)',1);
    yfit = P(1)*data_mi_input(:,10)'*69+P(2);
    hold on;
    plot(data_mi_input(:,10)'*69,yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w');
%% 
% cmap='plasma'
% fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 300])
% scatter(data_mi_input(:,10)*69,ba_input(:,6),30,tr,'filled');
% ylabel('Basal Branch Length_{mean}');xlabel('Input Horizontal extent L2/3_{EX} (µm)');set(gca,'FontSize',10);
%  hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
%  ref.YData=[90 0];;ref.Color='k';ylim([20 90]); ref.XData=[0 1200]
%  colormap(cmap);box off; hold on; c=colorbar;%c.Label.String = 'Relative pial position'
% %  text(-2,0.99,'r=0.44');text(-2,0.9,'p<0.001');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST FOR CLUSTERABILITY%%%%%%%%%%%%%%%







%% 



%% Passive ephys and apical
% ephys_sel=[];ephys_sel=[data_ephys(:,13:17)];
% morph_sel=[];morph_sel=[data_morpho(:,1:9) max_a' dis_peaka'];
% score_ephys=[];
% [coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
% [coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(zscore(morph_sel));
% var_exp(explained_ephys,[],[]); 
% var_exp(explained_morph,[],[]); 

%% Active ephys and basal
% ephys_sel_a=[];ephys_sel_a=[data_ephys(:,[1:7 18])];
% morph_sel_ba=[];morph_sel_ba=[data_morpho(:,10:19) max_b' dis_peakb'];
% [coeff_ephys_a,score_ephys_a,latent_ephys_a,~,explained_ephys_a,mu] = pca(zscore(ephys_sel_a));
% [coeff_morph_ba,score_morph_ba,latent_morph_ba,~,explained_morph_ba,mu] = pca(zscore(morph_sel_ba));
%% PC with dip test
dip_test_SW(score_morph_a(:,[1 2 3]),0,{'PC1','PC2','PC3'});
%% 

dip_test_SW(score_ephys(:,[1 2 3]),0,{'PC1','PC2','PC3'});
%% 

dip_test_SW(score_com(:,[1 2 3]),0,{'PC1','PC2','PC3'});
%% Ephys PC1 vs PC2
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250])
scatter(score_ephys(:,1),score_ephys(:,2),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('PC2 ephys');xlabel('PC1 ephys');set(gca,'FontSize',10);
%% Morpho PC1 vs PC2
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250])
scatter(score_morph(:,1),score_morph(:,2),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('PC2 morph');xlabel('PC1 morph');set(gca,'FontSize',10);
%% Input PC1 vs PC2
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250])
scatter(score_com(:,1),score_com(:,2),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('PC2 input');xlabel('PC1 input');set(gca,'FontSize',10);
%% 
tr=[]; tr=rescale(pia_ephys);
%umap_plot(ephys_sel, pia_ephys,'Pial depth');
umap_plot(score_ephys(:,[1 2 3]), tr,'');
%% 
umap_plot(score_ephys_a(:,[1 2 3]), pia_ephys,'Pial depth');
 
%% Mopho all properties apical
tr=[]; tr=rescale(pia_morpho);
%umap_plot([data_morpho(:,1:9) max_a' dis_peaka'], pia_morpho,'Pial depth')
umap_plot(score_morph_a(:,[1 2 3]), tr,'');
%% basal
umap_plot(score_morph_ba(:,[1 2 3]), pia_morpho,'Pial depth');
%% Input 
tr=[]; tr=rescale(pia_input');
umap_plot(score_com(:,[1 2 3]), tr,'Pial depth');



%% INVIVO

%% in vivo in vitro comparison
% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
directory=out_dir;% use cobined date structure named 
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 [od_out sftf_out sftf_sel sftf_pref spon_out pia_all delta_Ca fit_Ca] = concat_invivo(L23_PC_new);
 %%   %INVIVO INVITRO
od_out_iviv=[[str(imaps_id).OSIpref];[str(imaps_id).DSIpref];[str(imaps_id).ODIpref];[str(imaps_id).ORIpref];[str(imaps_id).DIRpref]...
              ;[str(imaps_id).Capeakpref];[str(imaps_id).Sigmapref];[str(imaps_id).SF];[str(imaps_id).TF];[str(imaps_id).sad];[str(imaps_id).noise];[str(imaps_id).pci]]';
          %% 
          
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
 
 
%% Sorted correlation checked with multiple comparison
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
%% 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 300])
imagesc(coeff_invivo([3 5 4 1 2],1:3));
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:3]);yticks([1:1:5]);colorbar;
xticklabels({'PC1','PC2','PC3'});
yticklabels({'ODI','R / R0_{max}','Tuning Width','gOSI','gDSI'});
set(gca,'FontSize',10);
%% Show PC2 vs pial depth 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 300])
scatter(score_invivo(:,2),invivo_feat(:,6),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 hold on;
 %ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[30 90];
% xticks([-5:2:5]);
%ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';
ylabel('Relative pial position');xlabel('PC2 (22%)')
%xlim([-5 5]);%xticks([40:80:280]);
text(4,0.99,'r=0.1');text(4,0.9,'p<0.01');
axis square
P=[];yfit=[];
 P = polyfit(score_invivo(:,2),invivo_feat(:,6),1);
    yfit = P(1)*score_invivo(:,2)'+P(2);
    hold on;
    plot(score_invivo(:,2)',yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 

%% PC1 invivo vs PC2 invivo 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250])
scatter(score_invivo(:,1),score_invivo(:,2),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('PC2 visresp');xlabel('PC1 visresp');set(gca,'FontSize',10);xticks([-5:2.5:5])
%% Dip test invivo using the first 3 PCs
dip_test_SW(score_invivo(:,[1 2 3]),0,{'PC1','PC2','PC3'});
%% UMAP plot
tr=rescale(pia_all);
umap_plot(score_invivo(:,[1 2 3]), invivo_feat(:,6)','');

%% 

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
scatter(od_out(:,3),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('Relative pial position');xlabel('ODI');set(gca,'FontSize',10);
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[0 1];
 yticks([0:0.25:1]);
 %ref.Color='k'; 
  %text(-2,0.99,'r=0.44');text(-2,0.9,'p<0.001');
  
   P=[];yfit=[];
   ty=find(isnan(od_out(:,3))==0);
 P = polyfit(od_out(ty,3),tr(ty)',1);
    yfit = P(1)*od_out(ty,3)+P(2);
    hold on;
    plot(od_out(ty,3),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
%% 

 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
scatter(max(delta_Ca(:,:),[],2),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('Relative pial position');xlabel('ODI');set(gca,'FontSize',10);
 %hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 %ref.YData=[0 1];
 yticks([0:0.25:1]);
 %ref.Color='k'; 
  %text(-2,0.99,'r=0.44');text(-2,0.9,'p<0.001');
  
   P=[];yfit=[];
   ty=find(isnan(od_out(:,2))==0);
 P = polyfit(od_out(ty,2),tr(ty)',1);
    yfit = P(1)*od_out(ty,2)+P(2);
    hold on;
    plot(od_out(ty,2),yfit,'-','Color',[0.5 0.5 0.5]);set(gca,'box','off');set(gcf,'color','w'); 
  %% Responiveness per depth in vivo all 
  tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.33);
g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.33 & tr'<=0.66);  
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.66); 

 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 350, 400]);
 subplot(3,1,1)
 h1=histogram(max(delta_Ca(g1,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.35]);h1.BinWidth = 20
hold on;p1=plot(nanmean(max(delta_Ca(g1,:),[],2)),0.35,'kv');p1.MarkerFaceColor='k'
hold on;p1=plot(nanmedian(max(delta_Ca(g1,:),[],2)),0.35,'mv');p1.MarkerFaceColor='m'
 title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(3,1,2)
 h1=histogram(max(delta_Ca(g2,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.35]);
  hold on;p1=plot(nanmean(max(delta_Ca(g2,:),[],2)),0.35,'kv');hold on;p1.MarkerFaceColor='k'
  p1=plot(nanmedian(max(delta_Ca(g2,:),[],2)),0.35,'mv');h1.BinWidth = 20;p1.MarkerFaceColor='m'
  title('>0.33 r. pia pos. <0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);ylabel('Percent of cells');set(gca,'FontSize',10)
  subplot(3,1,3)
 h1=histogram(max(delta_Ca(g3,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.35]);
 hold on;p1=plot(nanmean(max(delta_Ca(g3,:),[],2)),0.35,'kv');hold on;p1.MarkerFaceColor='k'
 p1=plot(nanmedian(max(delta_Ca(g3,:),[],2)),0.35,'mv');h1.BinWidth = 20;p1.MarkerFaceColor='m'
 title('r. pia pos. >0.66','FontWeight','Normal','FontSize',10);
xlabel('R / R0_{max}'); set(gca,'FontSize',10);
 
% g1=[];g2=[];g3=[];
% g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.25);
% g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'>200 & pia_all'<300);  
% g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.5);  
% par=[];
%   par=max(delta_Ca(:,:),[],2);
%   [statsout]=dual_barplot(par,g1,g3,2);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Overal binary responsiveness: no difference
tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(tr'<=0.5);
g2=find(tr'>0.25 & tr'<=0.5);  
g3=find(tr'>0.5); 

sum(od_out(g1,8))/length(od_out(g1,8))

sum(od_out(g2,8))/length(od_out(g2,8))

sum(od_out(g3,8))/length(od_out(g3,8))
%% Responiveness per depth in vivo all using sum
tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.33);
g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.33 & tr'<=0.66);  
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.66);  
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 400]);
subplot(3,1,1)
h1=histogram(sum(delta_Ca(g1,:),2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 600]);
ylim([0 0.3]);h1.BinWidth = 20;yticks([0:0.15:0.3]);set(gca,'FontSize',10)
hold on;p1=plot(nanmean(sum(delta_Ca(g1,:),2)),0.3,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(sum(delta_Ca(g1,:),2)),0.3,'mv');p1.MarkerFaceColor='m';
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(3,1,2)
h1=histogram(max(delta_Ca(g2,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 600]);
ylim([0 0.3]);ylabel('Percent of cells');set(gca,'FontSize',10);h1.BinWidth = 20;yticks([0:0.15:0.3]);set(gca,'FontSize',10)
hold on;p1=plot(nanmean(sum(delta_Ca(g2,:),2)),0.3,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(sum(delta_Ca(g2,:),2)),0.3,'mv');p1.MarkerFaceColor='m';
title('>0.33 r. pia pos. <0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
 subplot(3,1,3)
h1=histogram(sum(delta_Ca(g3,:),2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 600]);ylim([0 0.35]);
hold on;p1=plot(nanmean(sum(delta_Ca(g3,:),2)),0.3,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(sum(delta_Ca(g3,:),2)),0.3,'mv');p1.MarkerFaceColor='m';
ylim([0 0.3]);h1.BinWidth = 20;yticks([0:0.15:0.3]);set(gca,'FontSize',10);
title('r. pia pos. >0.66','FontWeight','Normal','FontSize',10);
xlabel('R / R0_{sum}'); set(gca,'FontSize',10);
%% Responiveness per depth in vivo all using sum, split L2/3 in half
g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.5);
g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'>200 & pia_all'<300);  
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.5);  
par=[];
%par=sum(delta_Ca(:,:),2);
par=max(delta_Ca(:,:),[],2);
[statsout]=dual_barplot(par,g1,g3,2);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
 
%% Violin boxplot;  split L2/3 in half: Responsiveness amplitude
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 260, 400]);
hold on;violin({par(g1) par(g3)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','m'); box off; ylabel('R / R0_{max}'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g3));
xlabel('Relative pial position');hold on;

%% 
figure;


hold on;pS=plotSpread({par(g1);par(g3)},'categoryIdx',[ones(length(par(g1)),1);ones(length(par(g3)),1)*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{'r','b'});hold on;

dat=[];
   dat=[par(g1);par(g3)];
   gro=[ones(length(par(g1)),1);ones(length(par(g3)),1)*2];
 
boxplot(dat,gro,'Colors','k','Symbol','');
  %% ODI per depth in vivo all 
tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.33);
g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.33 & tr'<=0.66);  
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.66);  
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 400]);
subplot(3,1,1)
h1=histogram(od_out(g1,3),20,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(od_out(g1,3)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(od_out(g1,3)),0.15,'mv');h1.BinWidth = 0.1;p1.MarkerFaceColor='m'
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);ylim([0 0.2]);
subplot(3,1,2);
h1=histogram(od_out(g2,3),20,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(od_out(g2,3)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(od_out(g2,3)),0.15,'mv');ylabel('Percent of cells');set(gca,'FontSize',10);p1.MarkerFaceColor='m'
h1.BinWidth = 0.1;ylim([0 0.2]);
title('>0.33 r. pia pos. <0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(3,1,3)
h1=histogram(od_out(g3,3),20,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(od_out(g3,3)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(od_out(g3,3)),0.15,'mv');xlabel('ODI');set(gca,'FontSize',10);p1.MarkerFaceColor='m';
h1.BinWidth =0.1;ylim([0 0.2]);
title('r. pia pos. >0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
%% stats barplot;  split L2/3 in half: ODI
g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'<=0.5);
g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'>200 & pia_all'<300);  
g3=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & tr'>0.5);  
par=[];
%par=abs(od_out(:,3));
par=(od_out(:,3));
[statsout]=dual_barplot(par,g1,g3,0);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Violin boxplot;  split L2/3 in half: Responsiveness amplitude
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 260, 400]);
par=[];
par=od_out(:,3);
hold on;violin({par(g1) par(g3)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','m'); box off; ylabel('ODI'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g3))
  xlabel('Relative pial position');
%% 
con_e= find(od_out(:,9)==1);
ips_e= find(od_out(:,10)==1);
bin_e= find(od_out(:,9)==1 & od_out(:,10)==1);

tr=[];tr=rescale(pia_all);g1c=[];g2c=[];g3c=[];
g1c=find(tr(con_e)'<=0.33);
g2c=find(tr(con_e)'>0.33 & tr(con_e)'<=0.66);  
g3c=find(tr(con_e)'>0.66); 
tr=[];tr=rescale(pia_all);g1i=[];g2i=[];g3i=[];
g1i=find(tr(ips_e)'<=0.33);
g2i=find(tr(ips_e)'>0.33 & tr(ips_e)'<=0.66);  
g3i=find(tr(ips_e)'>0.66); 
tr=[];tr=rescale(pia_all);g1b=[];g2b=[];g3b=[];
g1b=find(tr(bin_e)'<=0.33);
g2b=find(tr(bin_e)'>0.33 & tr(bin_e)'<=0.66);  
g3b=find(tr(bin_e)'>0.66); 

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 150, 400]);
yt=[length(g1c)/(length(g1c)+length(g1i)+length(g1b)); length(g1i)/(length(g1c)+length(g1i)+length(g1b)) ; length(g1b)/(length(g1c)+length(g1i)+length(g1b))]
subplot(3,1,1)
b1=bar(1,yt(1));;b1(1).FaceColor='b';hold on;b1=bar(2,yt(2));;b1(1).FaceColor='r';hold on;b1=bar(3,yt(3));;b1(1).FaceColor='w';box off;ylim([0 0.6]);set(gca,'FontSize',10);

subplot(3,1,2)
yt=[length(g2c)/(length(g2c)+length(g2i)+length(g2b)) length(g2i)/(length(g2c)+length(g2i)+length(g2b)) length(g2b)/(length(g2c)+length(g2i)+length(g2b))]
b1=bar(1,yt(1));;b1(1).FaceColor='b';hold on;b1=bar(2,yt(2));;b1(1).FaceColor='r';hold on;b1=bar(3,yt(3));;b1(1).FaceColor='w';box off; ylim([0 0.6]);ylabel('Fraction of cells');set(gca,'FontSize',10);

subplot(3,1,3)
yt=[length(g3c)/(length(g3c)+length(g3i)+length(g3b)) length(g3i)/(length(g3c)+length(g3i)+length(g3b)) length(g3b)/(length(g3c)+length(g3i)+length(g3b))]
b1=bar(1,yt(1));;b1(1).FaceColor='b';hold on;b1=bar(2,yt(2));;b1(1).FaceColor='r';hold on;b1=bar(3,yt(3));;b1(1).FaceColor='w';box off;ylim([0 0.6]);set(gca,'FontSize',10);
%% OSI
  tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & sum(delta_Ca,2)>0 & tr'<=0.5);
g2=find(od_out(:,8)==1 & sum(delta_Ca,2)>100 & pia_all'>200 & pia_all'<300);  
g3=find(od_out(:,8)==1 & sum(delta_Ca,2)>0 & tr'>0.5);  
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 300]);
subplot(2,1,1)
h1=histogram(od_out(g1,1),20,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;plot(nanmean(od_out(g1,1)),0.35,'kv');hold on;plot(nanmedian(od_out(g1,1)),0.35,'mv');
% subplot(3,1,2)
% h1=histogram(max(delta_Ca(g2,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.45]);
% hold on;plot(nanmean(max(delta_Ca(g2,:),[],2)),0.4,'kv');hold on;plot(nanmedian(max(delta_Ca(g2,:),[],2)),0.4,'mv');
 subplot(2,1,2)
h1=histogram(od_out(g3,1),20,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;plot(nanmean(od_out(g3,1)),0.35,'kv');hold on;plot(nanmedian(od_out(g3,1)),0.35,'mv');
%% 
%% stats barplot;  split L2/3 in half: OSI
qua = quantile(max(delta_Ca,[],2),0.75)
g1=[];g2=[];g3=[];tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & r2_all'>0.3 & max(delta_Ca,[],2)>qua & tr'<=0.5);
g2=find(od_out(:,8)==1 & r2_all'>0.3 & max(delta_Ca,[],2)>0 & pia_all'>200 & pia_all'<300);  
g3=find(od_out(:,8)==1 &  r2_all'>0.3 & max(delta_Ca,[],2)>qua & tr'>0.5);  
par=[];
par=od_out(:,1);
[statsout]=dual_barplot(par,g1,g3,0);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Violin boxplot;  split L2/3 in half: Responsiveness amplitude
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 260, 400]);
par=[];
par=od_out(:,1);
hold on;violin({par(g1) par(g3)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','m'); box off; ylabel('gOSI'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g3))
  xlabel('Relative pial position');

%  %% Responiveness per depth in vivo in vitro all
% tr=[];tr=rescale([str(imaps_id).pia_invivo]);g1=[];g2=[];g3=[];
% %tr=[];tr=rescale(pia_input);g1=[];g2=[];g3=[];
% % g1=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'<200);  
% % g2=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'>300);  
%  g1=find(od_out_iviv(:,6)>0 & tr'<0.5);  
%  g2=find(od_out_iviv(:,6)>0 & tr'>0.4);  
% fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 300]);
% subplot(2,1,1)
% h1=histogram(od_out_iviv(g1,6),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 350]);ylim([0 0.2]);
% hold on;plot(nanmean(od_out_iviv(g1,6)),0.2,'kv');hold on;plot(nanmedian(od_out_iviv(g1,6)),0.2,'mv');
% % subplot(3,1,2)
% % h1=histogram(max(delta_Ca(g2,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.45]);
% % hold on;plot(nanmean(max(delta_Ca(g2,:),[],2)),0.4,'kv');hold on;plot(nanmedian(max(delta_Ca(g2,:),[],2)),0.4,'mv');
%  subplot(2,1,2)
% h1=histogram(od_out_iviv(g2,6),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 350]);ylim([0 0.2]);
% hold on;plot(nanmean(od_out_iviv(g2,6)),0.2,'kv');hold on;plot(nanmedian(od_out_iviv(g2,6)),0.2,'mv');
%% 


%% Histogram for ORI PREF

g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & od_out(:,1)>0.25);
ORI_pref=od_out(g1,4);
% get the pia and also filter
piald = pia_all(g1);
piald=rescale(piald);

% ORI_pref(ORI_pref>180) = ORI_pref(ORI_pref>180) - 180;
ORI_pref = ORI_pref + 90;
ORI_pref(ORI_pref>179.999999) = ORI_pref(ORI_pref>179.999999) - 180;
angle_edges=[];
% separate the orientation in bins based on pia
[N,edges,bin] = histcounts(piald,[0 0.5  1]);
% [N,edges,bin] = histcounts(piald,3);
%angle_edges = [0 35 80 125 170 225]-22.5;
angle_edges = [0 45 90 135 180];
% allocate memory for the bin results
% depth_bins = zeros(size(N,2),size(angle_edges,2)-2);
depth_bins = zeros(size(N,2),size(angle_edges,2)-1);
% for all the depth bins
for bins = 1:size(N,2)
    % bin the angles
    temp_bins = histcounts(ORI_pref(bins==bin),angle_edges,'Normalization','probability');
    depth_bins(bins,:) = temp_bins
%     % combine the edge bins
%     depth_bins(bins,:) = [temp_bins(2:4),sum(temp_bins([1,5]))];
end

% plot the results
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 300]);
b1=bar(depth_bins');box off
%set(gca,'XTickLabels',string(angle_edges(2:5)+22.5))
set(gca,'XTickLabels',string(angle_edges));
legend({'<0.5','>0.5'});legend boxoff
xlabel('Preferred Orientation (deg)')
ylabel('Fraction of cells');set(gca,'FontSize',10);
b1(1).FaceColor='k';
b1(2).FaceColor=[0.5 0.5 0.5];
%% %% stats barplot;  split L2/3 in half: OSI
g1=[];g2=[];g3=[];tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & od_out(:,1)>0.25 & max(delta_Ca,[],2)>0 & tr'<=0.5);
g2=find(od_out(:,8)==1 & od_out(:,1)>0.25 & max(delta_Ca,[],2)>0 & pia_all'>200 & pia_all'<300);  
g3=find(od_out(:,8)==1 &  od_out(:,1)>0.25 & max(delta_Ca,[],2)>0 & tr'>0.5);  
par=[];
par=od_out(:,4);
[statsout]=dual_barplot(par,g1,g3,0);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Violin boxplot;  split L2/3 in half: Responsiveness amplitude
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 260, 400]);
par=[];
par=od_out(:,4);
hold on;violin({par(g1) par(g3)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','m'); box off; ylabel('Preffered Orientation (deg)'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g3))
  xlabel('Relative pial position');
  
%% 
ORI_pref=[];ORI_pref=od_out(:,4);
ORI_pref = ORI_pref + 90;
ORI_pref(ORI_pref>179.999999) = ORI_pref(ORI_pref>179.999999) - 180;

tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & od_out(:,1)>0.25 & tr'<=0.33);
g2=find(od_out(:,8)==1 & od_out(:,1)>0.25 & tr'>0.33 & tr'<=0.66);  
g3=find(od_out(:,8)==1 & od_out(:,1)>0.25 & tr'>0.66);  

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 400]);
subplot(3,1,1)
h1=histogram(ORI_pref(g1),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(ORI_pref(g1)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(ORI_pref(g1)),0.15,'mv');p1.MarkerFaceColor='m'
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(3,1,2);
h1=histogram(ORI_pref(g2),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(ORI_pref(g2)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(ORI_pref(g2)),0.15,'mv');p1.MarkerFaceColor='m'
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
title('>0.33 r. pia pos. <0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(3,1,3)
h1=histogram(ORI_pref(g3),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(ORI_pref(g3)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(ORI_pref(g3)),0.15,'mv');p1.MarkerFaceColor='m'
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
title('r. pia pos. >0.66','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
%% 
ORI_pref=[];ORI_pref=od_out(:,4);
ORI_pref = ORI_pref + 90;
ORI_pref(ORI_pref>179.999999) = ORI_pref(ORI_pref>179.999999) - 180;

tr=[];tr=rescale(pia_all);g1=[];g2=[];g3=[];
g1=find(od_out(:,8)==1 & od_out(:,1)>0.25 & tr'<=0.5);
g3=find(od_out(:,8)==1 & od_out(:,1)>0.25 & tr'>0.5); 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 400]);
subplot(2,1,1)
h1=histogram(ORI_pref(g1),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(ORI_pref(g1)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(ORI_pref(g1)),0.15,'mv');p1.MarkerFaceColor='m'
title('r. pia pos. <0.5','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
subplot(2,1,2);
h1=histogram(ORI_pref(g3),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 450]);ylim([0 0.35]);
hold on;p1=plot(nanmean(ORI_pref(g3)),0.15,'kv');p1.MarkerFaceColor='k';hold on;p1=plot(nanmedian(ORI_pref(g3)),0.15,'mv');p1.MarkerFaceColor='m'
title('r. pia pos. <0.33','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
title('r. pia pos. >0.5','FontWeight','Normal','FontSize',10);set(gca,'FontSize',10);
 %% Responiveness per depth in vivo in vitro all
tr=[];tr=rescale([str(imaps_id).pia_invivo]);g1=[];g2=[];g3=[];
tr=[];tr=rescale(pia_input');g1=[];g2=[];g3=[];
% g1=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'<200);  
% g2=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'>300);  
 g1=find(od_out_iviv(:,6)>0  & tr'<0.25);  
 g2=find(od_out_iviv(:,6)>0 &  tr'>0.65); 
 casum=[str(imaps_id).Casumpref];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 300]);
subplot(2,1,1)
h1=histogram(casum(g1),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 350]);ylim([0 0.2]);
hold on;plot(nanmean(casum(g1)),0.2,'kv');hold on;plot(nanmedian(casum(g1)),0.2,'mv');
% subplot(3,1,2)
% h1=histogram(max(delta_Ca(g2,:),[],2),20,'Normalization','probability');h1.FaceColor='w';box off;xlim([0 450]);ylim([0 0.45]);
% hold on;plot(nanmean(max(delta_Ca(g2,:),[],2)),0.4,'kv');hold on;plot(nanmedian(max(delta_Ca(g2,:),[],2)),0.4,'mv');
 subplot(2,1,2)
h1=histogram(casum(g2),8,'Normalization','probability');h1.FaceColor='w';box off;%xlim([0 350]);ylim([0 0.2]);
hold on;plot(nanmean(casum(g2)),0.2,'kv');hold on;plot(nanmedian(casum(g2)),0.2,'mv');

 %% 
g1=[];g2=[];
g1=find(od_out_iviv(:,6)>0 & tr'<0.25);  
g2=find(od_out_iviv(:,6)>0 & tr'>0.5); 
par=[];
  par=casum
  [statsout]=dual_barplot(par,g1,g2,2);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
%% Violin ODI plot in vivo in vitro all
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 400]);
par=[];
par=od_out_iviv(:,3);
hold on;violin({par(g1) par(g2)},'xlabel',{'< 0.5','> 0.5'},'facecolor',[1 1 1;1 1 1],'edgecolor','k',...
'mc','k',...
'medc','r--'); box off; ylabel('ODI'),legend('');legend boxoff;
[p k]=ranksum(par(g1),par(g2),'tail','left')
  %% 
  g1=[];g2=[];
g1=find(od_out_iviv(:,6)<10);  
g2=find(od_out_iviv(:,6)>180); 
  par=L4fr(:,1)
  [statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
  %% 
  
 figure;histogram(fit_Ca)
 tt=max(delta_Ca,[],2)

 %% 
 a=[];a=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & od_out(:,3)>0.5)
 %% 
 
  discretize_plot(pia_all(a),2,od_out(a,1),1)
  %% 
   discretize_plot(pia_all,3,od_out(:,8),1)
%% 
%% 
plot_avg_maps(str_imap,44,ex_map,in_map,pia_input,1,0,[]);

%% 
  g1=[];
  g2=[];
   a=[];
 a=find(od_out(:,8)==1 & abs(od_out(:,3))>0); 
  g1=find(pia_all(a)<200)  
  g2=find(pia_all(a)>290)
  par=od_out(a,4)
  [statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);hold on;set(gca,'FontSize',10);xtickangle(45);
  %% 
  a=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0); 
   g1=[];
  g2=[];
  g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'<200);  
  g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'>300);  
  
  
  figure;histogram(max(delta_Ca(g1,:),[],2),10,'Normalization','probability');hold on;histogram(max(delta_Ca(g2,:),[],2),10,'Normalization','probability')
  nanmean(max(delta_Ca(g1,:),[],2))
   nanmean(max(delta_Ca(g2,:),[],2))
   
  
   %% 
    g1=[];
  g2=[];
  g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>50 & pia_all'<200);  
  g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>50 & pia_all'>300);  
     figure;histogram(od_out(g1,1),8,'Normalization','probability');hold on;histogram(od_out(g2,1),8,'Normalization','probability')
  nanmean(od_out(g1,1))
   nanmean(od_out(g2,1))
   
   %% 
     g1=[];
  g2=[];
  g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'<200);  
  g2=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'>300);  
     figure;histogram(od_out(g1,4),8,'Normalization','probability');hold on;histogram(od_out(g2,4),8,'Normalization','probability')
  nanmean(od_out(g1,4))
   nanmean(od_out(g2,4))
   %% 
       g1=[];
  g2=[];
  g1=find(pia_all'<200);  
  g2=find(pia_all'>300);  
     figure;histogram(sftf_out(g1,2),8,'Normalization','probability');hold on;histogram(sftf_out(g2,2),8,'Normalization','probability')
  nanmean(od_out(g1,4))
   nanmean(od_out(g2,4))
   %% 
    a=find(od_out_iviv(:,1)>0)
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2)...
    span(a,:) pia_input(a) od_out_iviv(a,[1:end])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%% 
    a=find(od_out_iviv(:,1)>0)
com=[];com=[score_com(a,[1 2 3]) od_out_iviv(a,[1:end])]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%% 

figure;scatter(com(:,3),com(:,7))
   
   %% 
 
   %% 
   figure;scatter([str(imaps_id).Ca_peak_sftf],[str(imaps_id).pia_invivo])
   
   
   %% 

   
   %% 
   
   
   %% 
   figure;scatter([str(imaps_id).pia_invivo]',L4fr(:,1));
   
   figure;scatter(od_out_iviv(:,6),[str(imaps_id).pia_invivo]',20,L4fr(:,1),'filled');
   %% 
   
   g1=[];
  g2=[];
  g1=find(pia_input<200);  
  g2=find(pia_input>250);  
     figure;histogram(od_out_iviv(g1,4),8,'Normalization','probability');
     hold on;histogram(od_out_iviv(g2,4),8,'Normalization','probability')
  nanmean(od_out_iviv(g1,4))
   nanmean(od_out_iviv(g2,4))
   %% 
   
  figure;scatter(od_out_iviv(:,4),L4fr(:,1));
  [r p]=corrcoef(od_out_iviv(:,5),L4fr(:,1),'rows','pairwise')
  %% 
  
   a=find(od_out_iviv(:,2)>0);  
   g1=[];
  g2=[];
g1=find(pia_input<200)
g2=find(pia_input>220)
par=od_out_iviv(:,1)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
%% 
 a=find(od_out_iviv(:,1)>0);  
   g1=[];
  g2=[];
g1=find(pia_input<200)
g2=find(pia_input>250)
par=od_out_iviv(:,1)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);

%% 

figure;scatter(L4fr(:,1),od_out_iviv(:,1),20,pia_input,'filled')
%% 

   a=find(od_out(:,1)>0);  
   g1=[];
  g2=[];
g1=find(pia_all<200) ;
g2=find(pia_all>300) ;
par=od_out(:,1)
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);

%% 
dip_test_SW(data_morpho(:,[1 3 7 8 9 10 ]),1,{'PC1','PC2','PC3'});
%% 

clu_num = 3;
%pcs =[];
pcs     =[1 2 3 4 5 6];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input, clustering_input, leafOrder] = hca(data_morpho(:,[1 3 7 8 9 10]),0,'ward',clu_num,pia_morpho,1,0.75);%call function for clustering




%% 

