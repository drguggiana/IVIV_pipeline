%% L23 morpho/ephys structure
filename=uipickfiles('FilterSpec','C:\Users\simonw\Dropbox\Morph_ephys')%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% read out ephys only cells
data_ephys=[[mephys(:).Vrest]' [mephys(:).Rin]' [mephys(:).tau]' [mephys(:).Sag]' [mephys(:).Rheo]'...
               [mephys(:).Vmin]' [mephys(:).Vpeak]' [mephys(:).Vthresh]' [mephys(:).Vslope]' [mephys(:).Vhalf]'...
               [mephys(:).Vamp]' [mephys(:).AHP]' [mephys(:).APhwidth]'  [mephys(:).APfreq]'];
pia_ephys= [mephys(1:137).pia]' ;  
%% Pial depth distributions
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 200, 200]);
histogram(pia_ephys,6,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);
%% Correlation matrix ephys+pia
[R,P]=corrcoef([data_ephys pia_ephys],'rows','pairwise');
%% 

fig3= figure;set(fig3, 'Name', '');set(fig3, 'Position', [200, 300, 600, 120]);set(gcf,'color','w');
str={'RMP','Rin','Tau','Sag','Rheo','V_{min}','V_{peak}','V_{thresh}', 'Vslope_{max}','V_{half}','Spike_{amplitude}',...
    'AHP_{max}','Spike_{half width}','Max Freq'};

imagesc(R(15,1:14));
xticks([1:14])
xticklabels(str);
yticklabels({''});xtickangle(45);[cmap]=buildcmap('gm');colormap(cmap);colorbar;set(gca,'FontSize',10);
caxis([-0.4 0.4])
%% START
%Load structure with 147 cells used for the paper
str_invitro       = 'D:\Postdoc_Margrie\Projects\L23\structure\';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));