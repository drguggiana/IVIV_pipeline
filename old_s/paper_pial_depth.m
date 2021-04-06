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
if isnan(str_all(morph_id(i)).pialD) ==0    
pia_morpho(i)=str_all(morph_id(i)).pialD;
else
pia_morpho(i)=str_all(morph_id(i)).pia;    
end
end
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
%% Correlation with pial depth: APICAL using sorted correlation plot
stri={'Radial Dis._{max} (µm)','Total Length (µm)','Path Length_{max} (µm)','Branch Points','Branch Order_{max}','Branch Length_{mean} (µm)','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch (µm)'};
sorted_correlation([data_morpho(:,1:9) max_a' dis_peaka'],pia_morpho',stri,'k')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Apical');
%% Correlation with pial depth: BASAL using sorted correlation plot
stri={'Radial Dis._{max} (µm)','Total Length (µm)','Path Length_{max} (µm)','Branch Points','Branch Order_{max}','Branch Length_{mean} (µm)','Width / Height','Width','Height','Nr. Branches','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,10:19) max_b' dis_peakb'],pia_morpho',stri,'m')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Basal')
% %% Sliding winddow analysis APICAL
% [p_wha1 p_wha2] = sliding_window_comp(pia_morpho,data_morpho(:,3));
% %[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,max_a);
% %% Sliding winddow analysis BASAL
% [p_blb1 p_blb2] = sliding_window_comp(pia_morpho,data_morpho(:,15));
% [p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,data_morpho(:,13));
%% Plot the 1 corrleation each for Apical and BAsal
tr = rescale(pia_morpho);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 300, 500]);
subplot(2,1,1);
scatter(data_morpho(:,7)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 title('Apical');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Width / Height');
text(4,0.99,['r=-0.64']);text(4,0.9,['p<0.001']);set(gca,'FontSize',10);
subplot(2,1,2);
scatter(data_morpho(:,15),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 title('Basal');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[30 90];xticks([0:30:90]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Branch Length_{mean} (µm)')
xlim([25 85]);xticks([25:20:85]);ref.XData=[25 85];
text(27,0.99,['r=0.3']);text(27,0.9,['p<0.001']);set(gca,'FontSize',10);


%% Ephys



%% %% Read out which ones have ephys or not 
str_all=str;
for i=1:length(str_all)
    if isnan(str_all(i).Vmin)==0
    ephys(i)=1;    
    else
    ephys(i)=0;
    end
end
ephys_id=find(ephys==1);
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
sorted_correlation([data_ephys(:,13:17)],pia_ephys',stri,'k')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Passive');
%% Correlation with pial depth: Active using sorted correlation plot
stri={'APV_{min} (mV)','APV_{peak} (mV)','APV_{thresh} (mV)','APV_{slope} (mV)','APV_{half} (mV)','APV_{amp} (mV)','AHP (mV)','APfreq (Hz)'}
sorted_correlation([data_ephys(:,[1:7 18])],pia_ephys',stri,[0 0.5 0.5])
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Active');
% %% Sliding winddow analysis Passive
% [p_wha1 p_wha2] = sliding_window_comp(pia_ephys,data_ephys(:,17));
% %% Sliding winddow analysis Active
% [p_blb1 p_blb2] = sliding_window_comp(pia_ephys,data_ephys(:,4));
%% Plot the 1 corrleation each for passive and active
tr=[];tr = rescale(pia_ephys);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_ephys(:,14)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor',[0 0.5 0.5]);
 %set(gca, 'YDir','reverse');
 title('Passive');text(55,0.99,'r=-0.29');text(55,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Tau (ms)')
xlim([10 70]);ref.XData=[15 70];xticks([10:20:80]);
subplot(1,2,2);
scatter(data_ephys(:,4),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
  title('Active');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[30 90];xticks([0:30:90]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('APV_{slope} (mV)')
xlim([40 280]);ref.XData=[40 280];xticks([40:80:280]);
text(240,0.19,'r=-0.31');text(240,0.1,'p<0.001');
%% MORPHO AND EPHYS
%% Overlapping morpho ephys
[~,ia,ib] =intersect(ephys_id,morph_id)
%% Read out 33 cells that overlap
data_morph_both=[];
data_morph_both=[data_morpho(ib,:)];
data_ephys_both=[data_ephys(ia,:)];
pia_both=[pia_morpho(ib)];
%% Correlation of both
%Passive with apical parameters
ap_morph=[data_morph_both(:,1:9) max_a(ib)' dis_peaka(ib)'];
%G2=correlation_matrix([ap_morph data_ephys_both(:,[13:17])],0);

[R2,P2]=corrcoef([ap_morph data_ephys_both(:,[13:17])],'rows','pairwise');

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 450, 300])
imagesc(R2(12:end,1:11));c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:11]);yticks([1:1:5])
yticklabels({'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'});set(gca,'FontSize',10);
xticklabels({'Radial Dis._{max}','Total Length','Path Length_{max}','Branch Points','Branch Order_{max}','Branch Length_{mean}','Width / Height','Width','Height','Peak Nr. cross','Dis. peak branch'});
xtickangle(45);set(gca,'FontSize',10);


%% PCA for morpho ephys

ephys_sel=data_ephys(:,[13 14 15 17]);
morph_sel=[data_morpho(:,[1 3 5 6 7 8 9]) max_a' dis_peaka'];
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
[coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(zscore(morph_sel));
var_exp(explained_ephys,[],[]); 
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
% cell
plot_avg_maps(str_imap,110,ex_map,in_map,pia_input,1,0,[]);
%% 
plot_avg_maps(str_imap,108,ex_map,in_map,pia_input,1,0,[]);
%% 
plot_avg_maps(str_imap,80,ex_map,in_map,pia_input,1,0,[]);
%% %Pial depth correlation with EX IN L23, L4, L5 for panel F
plot_fraction(L23fr,L4fr,L5fr,pia_input);
%% 
%% Correlation with pial depth: Active using sorted correlation plot
stri={'L2/3 Vfr_{EX}','L4 Vfr_{EX}','L5 Vfr_{EX}','L2/3 Vfr_{IN}','L4 Vfr_{IN}','L5 Vfr_{IN}','L2/3 Vfr_{d}','L4 Vfr_{d}','L5 Vfr_{d}',...
   'L2/3 HEX_{EX}','L4 HEx_{EX}','L5 HEx_{EX}','L2/3 HEx_{IN}','L4 HEx_{IN}','L5 HEx_{IN}','L2/3 HEx_{d}','L4 HEx_{d}','L5 HEx_{d}'}
sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2)...
    diffL23fr,diffL4fr,diffL4fr,span,spandL23,spandL4,spandL5],pia_input,stri,[0.5 0.5 0.5])
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Active');

%% 
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(L4fr(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');
 %set(gca, 'YDir','reverse');
 xlim([0 0.7])
 text(0.1,0.99,'r=-0.41');text(0.1,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Fraction')
ref.XData=[0 0.7];xticks([0:0.15:0.7]);

subplot(1,3,2);
scatter(L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
 %set(gca, 'YDir','reverse');
 xlim([0 0.6])
 text(0.1,0.99,'r=-0.41');text(0.1,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';;xlabel('Fraction')
ref.XData=[0 0.6];xticks([0:0.15:0.6]);

subplot(1,3,3);
scatter(L4fr(:,1)'-L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 xlim([-0.2 0.6]);
 xlabel('\Delta EX-IN vertical fraction')
%  hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
% ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Fraction')
% ref.XData=[0 0.6];xticks([0:0.15:0.6]);

%% 
sort_fractions_input([pia_input L4fr],1,'L2/3')



%% 
%% 
plot_avg_maps(str_imap,44,ex_map,in_map,pia_input,1,0,[]);

%% 
dip_test_SW(data_ephys(:,[4 13:15 17]),1,{'PC1','PC2','PC3'});
%% 
dip_test_SW(data_morpho(:,[1 3 7 8 9 10 ]),1,{'PC1','PC2','PC3'});
%% 

clu_num = 3;
%pcs =[];
pcs     =[1 2 3 4 5 6];
%pcs     =[5];
%including the 6 PCs = pial depth 
[idx_input, clustering_input, leafOrder] = hca(data_morpho(:,[1 3 7 8 9 10]),0,'ward',clu_num,pia_morpho,1,0.75);%call function for clustering