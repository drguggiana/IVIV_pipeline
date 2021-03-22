%% L23 morpho/ephys structure
filename=uipickfiles('FilterSpec','C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Stage4_allcells_structure')%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% Read out which ones have morphology or not 
str_all=str_allcells;
for i=1:length(str_all)
    if isempty(str_all(i).RDA)==0
    morph(i)=1;    
    else
    morph(i)=0;
    end
end
morph_id=find(morph==1);
%% Read out morphology parameters for apical and basal tree
  data_morpho=[[str_all(:).RDA]' [str_all(:).LA]' [str_all(:).PLA]' [[str_all(:).BPA]'...
               [str_all(:).BOA]' [str_all(:).BLA]' [str_all(:).WHA]' [str_all(:).XSA]' [str_all(:).YSA]' ...
               [str_all(:).RDB]' [str_all(:).LB]' [str_all(:).PLB]' [str_all(:).BPB]'...
               [str_all(:).BOB]' [str_all(:).BLB]' [str_all(:).WHB]' [str_all(:).XSB]' [str_all(:).YSB]' [str_all(:).NB]']];
%% Read out pial depth
for i=1:length(morph_id)
if isempty(str_all(morph_id(i)).pialD) ==0    
pia_morpho(i)=str_all(morph_id(i)).pialD;
else
pia_morpho(i)=str_all(morph_id(i)).pia;    
end
end
%% Sholl Analysis APICAL
for i=1:length(morph_id)
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(str_all(morph_id(i)).tr_apical, 20, '-s');
max_a(i)=max(s);
temp=dd(find(s==max(s)));
dis_peaka(i)=temp(1);
end
%% Sholl Analysis BASAL
for i=1:length(morph_id)
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(str_all(morph_id(i)).tr_basal, 20, '-s');
max_b(i)=max(s);
temp=dd(find(s==max(s)));
dis_peakb(i)=temp(1);
end
%% Correlation with pial depth: APICAL using sorted correlation plot
stri={'RD_{max} (µm)','Total Length (µm)','PL_{max} (µm)','BP','BO','BL (µm)','Width/Height','Width','Height','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,1:9) max_a' dis_peaka'],pia_morpho',stri,'m')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Apical')
%% Correlation with pial depth: BASAL using sorted correlation plot
stri={'RD_{max} (µm)','Total Length (µm)','PL_{max} (µm)','BP','BO','BL (µm)','Width/Height','Width','Height','BTN','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,10:19) max_b' dis_peakb'],pia_morpho',stri,'k')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Basal')
%% Sliding winddow analysis APICAL
[p_wha1 p_wha2] = sliding_window_comp(pia_morpho,data_morpho(:,3));
%[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,max_a);
%% Sliding winddow analysis BASAL
[p_blb1 p_blb2] = sliding_window_comp(pia_morpho,data_morpho(:,15));
[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,data_morpho(:,13));
%% 

[cda] = RAFisher2cda(data_morpho(:,2)',1,1,0.05)
%% Plot the 1 corrleation each for Apical and BAsal
tr = rescale(pia_morpho);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_morpho(:,7)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Width / Height')
subplot(1,2,2);
scatter(data_morpho(:,15),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[30 90];xticks([0:30:90]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Mean Branch length (µm)')
xlim([30 85]);ref.XData=[30 85];

%% Ephys
%% %% Read out which ones have ephys or not 
str_all=str_allcells;
for i=1:length(str_all)
    if isempty(str_all(i).Vmin)==0
    ephys(i)=1;    
    else
    ephys(i)=0;
    end
end
ephys_id=find(ephys==1);
%% Read out morphology parameters for apical and basal tree
  data_ephys=[[str_all(:).Vmin]' [str_all(:).Vpeak]' [str_all(:).Vthresh]' [[str_all(:).Vslope]'...
               [str_all(:).Vhalf]' [str_all(:).Vamp]' [str_all(:).AHP]' [str_all(:).APrise]' [str_all(:).APfall]' ...
               [str_all(:).APbwidth]' [str_all(:).APhwidth]' [str_all(:).APlate]' [str_all(:).Vrest]'...
               [str_all(:).tau]' [str_all(:).Rin]' [str_all(:).Sag]' [str_all(:).Rheo]' [str_all(:).APfreq]']];

%% %% Read out pial depth
for i=1:length(ephys_id)
pia_ephys(i)=str_all(ephys_id(i)).pia;    
end

%% Correlation with pial depth: Passive using sorted correlation plot
stri={'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Sag Ratio','Rheobase'}
sorted_correlation([data_ephys(:,13:17)],pia_ephys',stri,'g')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Passive');

%% Correlation with pial depth: Active using sorted correlation plot
stri={'V_{min} (mV)','V_{peak} (mV)','V_{thresh} (mV)','V_{slope} (mV)','V_{half} (mV)','V_{amp} (mV)','AHP (mV)','APfreq (Hz)'}
sorted_correlation([data_ephys(:,[1:7 18])],pia_ephys',stri,'k')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);hold on;title('Active AP');
%% Sliding winddow analysis Passive
[p_wha1 p_wha2] = sliding_window_comp(pia_ephys,data_ephys(:,17));
%% Sliding winddow analysis Active
[p_blb1 p_blb2] = sliding_window_comp(pia_ephys,data_ephys(:,4));
%% Plot the 1 corrleation each for passive and active
tr=[];tr = rescale(pia_ephys);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_ephys(:,14)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','g');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Tau (ms)')
xlim([10 75]);ref.XData=[15 75];
subplot(1,2,2);
scatter(data_ephys(:,4),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; %ref.XData=[30 90];xticks([0:30:90]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('AP slope')
xlim([50 250]);ref.XData=[50 250];
%% 
%% Overlapping morpho ephys
[~,ia,ib] =intersect(ephys_id,morph_id)
%% REad out 33 cells that overlap
data_morph_both=[data_morpho(ib,:)];
data_ephys_both=[data_ephys(ia,:)];
pia_both=[pia_morpho(ib)];
%% 
G2=correlation_matrix([data_morph_both data_ephys_both(:,[1:7 18])],0);

G2=correlation_matrix([data_morph_both data_ephys_both(:,[13:17])],0);
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