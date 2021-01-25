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
%% Sholl Analysis


%% Correlation with pial depth: APICAL
sorted_correlation(data_morpho(:,1:9),pia_morpho')


%% 

hold on;
tidx=[];
tidx=find(P1(end,9:end-1)<0.05);
corr_basal=R1(end,10:end-1);
hold on;plot(10:19,corr_basal,'k.')
hold on;plot(tidx+9,corr_basal(tidx),'k*');

%hold on;line([1 size(data_morpho,2)], [0.05 0.05],'Color','k','LineStyle','--');
%% 
sliding_window_comp(pia_morpho,data_morpho(:,3))
%% 
       