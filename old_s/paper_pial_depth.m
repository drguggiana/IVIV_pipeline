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
scatter(data_ephys(:,14)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 title('Passive');text(55,0.99,'r=-0.29');text(55,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Tau (ms)')
xlim([10 70]);ref.XData=[15 70];xticks([10:20:80]);
subplot(1,2,2);
scatter(data_ephys(:,4),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor',[0 0.5 0.5]);
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
var_exp(explained_morph,[],[]); 
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
%% Correlation with pial depth: Vertical fraction
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
sorted_correlation([L23fr(:,1),L4fr(:,1),L5fr(:,1),L23fr(:,2),L4fr(:,2),L5fr(:,2),diffL23fr,diffL4fr,diffL4fr],pia_input,stri,[0.5 0.5 0.5])
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Vertical fraction');
%% Correlation with pial depth: Maximum horizontal extent
stri={'L2/3_{EX}','L4_{EX}','L5_{EX}','L2/3_{IN}','L4_{IN}','L5_{IN}','L2/3_{EX-IN}','L4_{EX-IN}','L5_{EX-IN}'}
sorted_correlation([span,spandL23,spandL4,spandL5],pia_input,stri,[0.5 0.5 0.5])
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Maximal horizontal extent');
%% Plot most significant corrleations vertical fractions
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(L4fr(:,1)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');
 %set(gca, 'YDir','reverse');
 xlim([0 0.7])
 text(0.1,0.99,'r=0.41');text(0.1,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Fraction')
ref.XData=[0 0.7];xticks([0:0.15:0.7]);

subplot(1,3,2);
scatter(L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
 %set(gca, 'YDir','reverse');
 xlim([0 0.6])
 text(0.1,0.99,'r=0.3');text(0.1,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';;xlabel('Fraction')
ref.XData=[0 0.6];xticks([0:0.15:0.6]);

subplot(1,3,3);
scatter(L4fr(:,1)'-L4fr(:,2)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 xlim([-0.2 0.6]);
 xlabel('\Delta EX-IN vertical fraction');
 hold on;line([0 0], [0 1],'Color','k','LineStyle','--');xlim([-0.4 0.5]);
%  hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
% ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Fraction')
% ref.XData=[0 0.6];xticks([0:0.15:0.6]);

%% Plot most significant corrleations horizontal extent
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(span(:,1)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');
 %set(gca, 'YDir','reverse');
 xlim([0 15*69]);
 text(800,0.99,'r=-0.3');text(800,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Horizontal Extent (µm)')
 ref.XData=[0 15*69];xticks([0:500:1000]);

subplot(1,3,2);
scatter(span(:,4)'*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
 %set(gca, 'YDir','reverse');
 xlim([0 15*69]);
 text(800,0.99,'r=-0.22');text(800,0.9,'p<0.001');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Horizontal Extent (µm)')
 ref.XData=[0 15*69];xticks([0:500:1000]);

subplot(1,3,3);
scatter(spandL23*69,tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 xlim([-500 500]);
 xlabel('\Delta EX-IN vertical fraction');

 %ref.XData=[0 4];xticks([0:1:4]);
ylabel('Relative pial position');xlabel('\Delta EX-IN horizontal extent')
 hold on;line([0 0], [0 1],'Color','k','LineStyle','--');set(gca,'FontSize',10)
%% PCA with whole maps
 [coeff_ex,score_ex,coeff_in,score_in,coeff_com,score_com] = map_align_PCA(str,imaps_id);
%% Plot correlation matrix of map scores
tr=[];tr = rescale(pia_input);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
scatter(score_com(:,1),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
ylabel('Relative pial position');xlabel('PC1_{com}');set(gca,'FontSize',10);
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
 ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k'; 
  text(-2,0.99,'r=0.44');text(-2,0.9,'p<0.001');
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST FOR CLUSTERABILITY%%%%%%%%%%%%%%%







%% 



%% Passive ephys and apical
ephys_sel=[];ephys_sel=[data_ephys(:,13:17)];
morph_sel=[];morph_sel=[data_morpho(:,1:9) max_a' dis_peaka'];
score_ephys=[];
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca(zscore(ephys_sel));
[coeff_morph,score_morph,latent_morph,~,explained_morph,mu] = pca(zscore(morph_sel));
var_exp(explained_ephys,[],[]); 
var_exp(explained_morph,[],[]); 

%% Active ephys and basal
ephys_sel_a=[];ephys_sel_a=[data_ephys(:,[1:7 18])];
morph_sel_ba=[];morph_sel_ba=[data_morpho(:,10:19) max_b' dis_peakb'];
[coeff_ephys_a,score_ephys_a,latent_ephys_a,~,explained_ephys_a,mu] = pca(zscore(ephys_sel_a));
[coeff_morph_ba,score_morph_ba,latent_morph_ba,~,explained_morph_ba,mu] = pca(zscore(morph_sel_ba));
%% PC with dip test
dip_test_SW(score_ephys(:,[1 2 3]),0,{'PC1','PC2','PC3'});
dip_test_SW(score_morph(:,[1 2 3]),0,{'PC1','PC2','PC3'});
dip_test_SW(score_com(:,[1 2 3]),0,{'PC1','PC2','PC3'});
%% 

%umap_plot(ephys_sel, pia_ephys,'Pial depth');
umap_plot(score_ephys(:,[1 2 3]), pia_ephys,'Pial depth');
%% 
umap_plot(score_ephys_a(:,[1 2 3]), pia_ephys,'Pial depth');
 
%% Mopho all properties apical
%umap_plot([data_morpho(:,1:9) max_a' dis_peaka'], pia_morpho,'Pial depth')
umap_plot(score_morph(:,[1 2 3]), pia_morpho,'Pial depth');
%% basal
umap_plot(score_morph_ba(:,[1 2 3]), pia_morpho,'Pial depth');
%% Input 
umap_plot(score_com(:,[1 2 3]), pia_input','Pial depth');



%% INVIVO

%% in vivo in vitro comparison
% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% %% Read out ori pref of all 
 [od_out sftf_out sftf_sel sftf_pref spon_out pia_all delta_Ca fit_Ca] = concat_invivo(L23_PC_new);
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
   
   %% Responiveness per depth
   
    a=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0); 
   g1=[];
  g2=[];
  g1=find(od_out(:,8)==1 & max(delta_Ca,[],2)>0 & pia_all'<200);  
  
 figure;histogram(max(delta_Ca(g1,:),[],2),10,'Normalization','probability');
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250]);
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
   %INVIVO INVITRO
   od_out_iviv=[[str(imaps_id).OSIpref];[str(imaps_id).DSIpref];[str(imaps_id).ODIpref];[str(imaps_id).ORIpref];[str(imaps_id).DIRpref]...
              ;[str(imaps_id).Capeakpref];[str(imaps_id).Sigmapref];[str(imaps_id).SF];[str(imaps_id).TF];[str(imaps_id).sad];[str(imaps_id).noise];[str(imaps_id).pci]]';
   
   %% 
   figure;scatter([str(imaps_id).Ca_peak_sftf],[str(imaps_id).pia_invivo])
   
   
   %% 
    g1=[];
  g2=[];
  g1=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'<200);  
  g2=find(od_out_iviv(:,6)>0 & [str(imaps_id).pia_invivo]'>250 & [str(imaps_id).pia_invivo]'<300);  
     figure;histogram(od_out_iviv(g1,6),8,'Normalization','probability');hold on;histogram(od_out_iviv(g2,6),8,'Normalization','probability')
  nanmean(od_out_iviv(g1,6))
   nanmean(od_out_iviv(g2,6))
   
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