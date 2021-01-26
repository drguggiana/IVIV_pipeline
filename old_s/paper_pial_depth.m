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
 stri={'RD_{max} (µm)','Total Length (µm)','PL_{max} (µm)','BP','BO','BL (µm)','WH','XS','YS','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,1:9) max_a' dis_peaka'],pia_morpho',stri,'m')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10)
%% Correlation with pial depth: BASAL using sorted correlation plot
stri={'RD_{max} (µm)','Total Length (µm)','PL_{max} (µm)','BP','BO','BL (µm)','WH','XS','YS','BTN','Peak Nr. cross','Dis. peak branch (µm)'}
sorted_correlation([data_morpho(:,10:19) max_b' dis_peakb'],pia_morpho',stri,'k')
hold on;xlabel('Correlation with pial depth')
hold on;set(gca,'FontSize',10);
%% Sliding winddow analysis APICAL
[p_wha1 p_wha2] = sliding_window_comp(pia_morpho,data_morpho(:,8));
%[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,max_a);
%% Sliding winddow analysis BASAL
[p_blb1 p_blb2] = sliding_window_comp(pia_morpho,data_morpho(:,15));
[p_max_a1 p_max_a2] = sliding_window_comp(pia_morpho,data_morpho(:,13));
%% 

%% 
tr = rescale(pia_morpho);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 200]);
subplot(1,2,1);
scatter(data_morpho(:,7)',tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','m');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[0 4];xticks([0:1:4]);
ref.YData=[1 0];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Width/Height')
subplot(1,2,2);
scatter(data_morpho(:,15),tr,20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 %set(gca, 'YDir','reverse');
 hold on;ref= refline(0,1);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[30 90];xticks([0:30:90]);
ref.YData=[0 1];yticks([0:0.25:1]);ref.Color='k';ylabel('Relative pial position');xlabel('Branch length')

%% 



%% 
figure;scatter(data_morpho(:,15),pia_morpho);
figure;scatter(data_morpho(:,13),pia_morpho);


%% 
       