%% L23 morpho/ephys structure
filename=uipickfiles('FilterSpec','C:\Users\simonw\Dropbox\Morph_ephys')%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% read out apical and basal traces for sholl analysis
close all;
for i=1:length(mephys)
    if ~isempty(mephys(i).RDA)==1 
    zz{:,i}=mephys(i).tr_apical;
    zz_b{:,i}=mephys(i).tr_basal;
    else
        zz{:,i}=NaN;
        zz_b{:,i}=NaN;
    end
end
%% sholl analysis
zz_both=[zz; zz_b];
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis2(zz_both,1:293,1);

%% Morphology
for i=1:length(mephys)
    if isempty(mephys(i).RDA)==0
    morph(i)=1;
    else
    morph(i)=0;
    end
end
  data_morpho=[[mephys(:).RDA]' [mephys(:).LA]' [mephys(:).PLA]' [[mephys(:).BPA]'...
               [mephys(:).BOA]' [mephys(:).BLA]' [mephys(:).WHA]' [mephys(:).XSA]' [mephys(:).YSA]' [max_s(~isnan(max_s))]' [dis_s(~isnan(max_s))]'  ...
               [mephys(:).RDB]' [mephys(:).LB]' [mephys(:).PLB]' [mephys(:).BPB]'...
               [mephys(:).BOB]' [mephys(:).BLB]' [mephys(:).WHB]' [mephys(:).XSB]' [mephys(:).YSB]' [mephys(:).NB]'] [max_s_ba(~isnan(max_s))]' [dis_s_ba(~isnan(max_s))]'];
  pia_morpho= [mephys(find(morph==1)).pia]';
%% Correlation with pia
G2=correlation_matrix([data_morpho pia_morpho],0);
%% Pltting significant Rs with pia
[R,P]=corrcoef([data_morpho pia_morpho],'rows','pairwise');
fig3= figure;set(fig3, 'Name', '');set(fig3, 'Position', [200, 300, 600, 120]);set(gcf,'color','w');
st={'RDA','LA','PLA','BPA','BOA','BLA','WHA','XSA','YSA','PSCA','PSDA', 'RDB','LB','PLB',...
    'BPB','BOB','BLB','WHB','XSB','YSB','NB','PSCB','PSDB'};
imagesc(G2(24,:));xticks([1:23]);xticklabels(st);
yticklabels({''});xtickangle(45);[cmap]=buildcmap('gwm');colormap(cmap);colorbar;set(gca,'FontSize',10);caxis([-1 1]);
%% PCA morpho 
[coeff_mo,score_mo,latent_mo,~,explained_mo,mu] = pca([zscore(data_morpho(:,:))]);
G2=correlation_matrix([score_mo(:,:) pia_morpho],0);
%% Display variance explained for ex and in maps
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y2=explained_mo;yyaxis left;h2=bar(y2);h2.EdgeColor = 'k';h2.FaceColor = [0 0 0.7];%ylim([0 20])
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:10:100]);
y4=cumsum(y2);yyaxis right
hold on;p2=plot(y4,'--')
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[1 0 0];p2.Color=[0 0 1];
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]);
%% Plot weights
figure;bar(coeff_mo(:,2),'k');set(gcf,'color','w');box off;ylabel('PC1 weight');xlabel('Feature number');set(gca,'FontSize',14);
%% PCA correlation with pial depth
corr_plot(score_mo(:,1),pia_morpho,[],{'PC1','Pial depth'})
%% Dip test
dip_test_SW(data_morpho,1,{'PC1','PC2','PC3'});
%% read out ephys only cells
data_ephys=[[mephys(:).Vrest]' [mephys(:).Rin]' [mephys(:).tau]' [mephys(:).Sag]' [mephys(:).Rheo]'...
               [mephys(:).Vmin]' [mephys(:).Vpeak]' [mephys(:).Vthresh]' [mephys(:).Vslope]' [mephys(:).Vhalf]'...
               [mephys(:).Vamp]' [mephys(:).AHP]' [mephys(:).APhwidth]'  [mephys(:).APfreq]'];
pia_ephys= [mephys(1:137).pia]' ;  
%% Pial depth distributions
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 200, 200]);
histogram(pia_ephys,6,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 70]);
fig2=figure;set(gcf,'color','w');set(fig2, 'Position', [100, 600, 200, 200]);
histogram(pia_morpho,6,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 70]);
%% Correlation matrix ephys+pia
[R,P]=corrcoef([data_ephys pia_ephys],'rows','pairwise');
G3=correlation_matrix([data_ephys(:,:) pia_ephys],0);
%% Correlation with pial depth
fig3= figure;set(fig3, 'Name', '');set(fig3, 'Position', [200, 300, 600, 200]);set(gcf,'color','w');
st={'RMP','Rin','Tau','Sag','Rheo','V_{min}','V_{peak}','V_{thresh}', 'Vslope_{max}','V_{half}','Spike_{amplitude}',...
    'AHP_{max}','Spike_{half width}','Max Freq'};
imagesc(G3(15,:));xticks([1:14]);xticklabels(st);
yticklabels({''});xtickangle(45);[cmap]=buildcmap('gwm');colormap(cmap);colorbar;set(gca,'FontSize',10);caxis([-1 1]);
%% PCA ephys passive
[coeff_ephys,score_ephys,latent_ephys,~,explained_ephys,mu] = pca([zscore(data_ephys(:,:))]);
G2=correlation_matrix([score_ephys(:,:) pia_ephys],0);
%% Dip test
dip_test_SW(data_ephys,1,{'PC1','PC2','PC3'});
%% 

%% Overlapping morpho ephys
eph=1:137;
morph_id=find(morph==1);
[~,ia,ib] =intersect(eph,find(morph==1));

%% Data for both
data_both=[[mephys(ia).Vrest]' [mephys(ia).Rin]' [mephys(ia).tau]' [mephys(ia).Sag]' [mephys(ia).Rheo]'...
               [mephys(ia).Vmin]' [mephys(ia).Vpeak]' [mephys(ia).Vthresh]' [mephys(ia).Vslope]' [mephys(ia).Vhalf]'...
               [mephys(ia).Vamp]' [mephys(ia).AHP]' [mephys(ia).APhwidth]'  [mephys(ia).APfreq]'...
               [mephys(ia).RDA]' [mephys(ia).LA]' [mephys(ia).PLA]' [[mephys(ia).BPA]'...
               [mephys(ia).BOA]' [mephys(ia).BLA]' [mephys(ia).WHA]' [mephys(ia).XSA]' [mephys(ia).YSA]' max_s(ia)' dis_s(ia)'...
               [mephys(ia).RDB]' [mephys(ia).LB]' [mephys(ia).PLB]' [mephys(ia).BPB]'...
               [mephys(ia).BOB]' [mephys(ia).BLB]' [mephys(ia).WHB]' [mephys(ia).XSB]' [mephys(ia).YSB]' [mephys(ia).NB]' max_s_ba(ia)' dis_s_ba(ia)']];
           
           
 %% Corrleation all
G4=correlation_matrix([data_both],0);

%% Tau vs apical length
corr_plot(data_both(:,3),data_both(:,16),[],{'Tau','Total length apical'})
%% Tau vs branch points apical
corr_plot(data_both(:,3),data_both(:,18),[],{'Tau','Branch points apical'})
%% AHP vs apical length
corr_plot(data_both(:,12),data_both(:,16),[],{'AHP','Total length apical'})
%% 
corr_plot(data_both(:,10),data_both(:,23),[],{'Vhalf','WHA'})
%% Sag vs basal tree number
corr_plot(data_both(:,4),data_both(:,35),[],{'Sag','Number basl trees'})
%% Sag vs basal tree number
corr_plot(data_both(:,4),data_both(:,29),[],{'Sag','Number basl trees'})
%% Rinput and YSA
corr_plot(data_both(:,2),data_both(:,23),[],{'Input resistance','Apical max vertical extent '});set(gca,'FontSize',10)
%% 
corr_plot(data_both(:,1),data_both(:,21),[],{'Vrest (mV)','Width / Height Apical'});set(gca,'FontSize',10)
%% AP max frequency vs apucal max peak crossing distance
corr_plot(data_both(:,14),data_both(:,25),[],{'Spike max frequency','Dis peak branch'})
%% Input maps

%% START
%Load structure with 147 cells used for the paper
str_invitro       = 'D:\Postdoc_Margrie\Projects\L23\structure\';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% 
%% %% Read out important parameters
% get a vector with the cells to use
iviv_cells = find([str(:).iviv]==1);
%Get normalized maps + raw maps
for i=1:length(str)
ex_map(:,:,i) = str(i).subpixel_excMap;
in_map(:,:,i) = str(i).subpixel_inhMap;
end
for i=1:length(str)
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
% od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
%              ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci];[str(:).error_pref];[str(:).r2_pref];[str(:).tun_pref]]';
od_out_iviv=[[str(:).OSIpref];[str(:).DSIpref];[str(:).ODIpref];[str(:).ORIpref];[str(:).DIRpref]...
              ;[str(:).Capeakpref];[str(:).Sigmapref];[str(:).SF];[str(:).TF];[str(:).sad];[str(:).noise];[str(:).pci]]';
%ipsi, contra, bino, unres, resp
 ipsi_id=find([str(:).ipsi]==1);
 contra_id=find([str(:).contra]==1);
 bino_id=find([str(:).bino]==1);
 unres_id=find([str(:).unres]==1);
 resp_id=find([str(:).resp]==1);
 %% 
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 200, 200]);
histogram(pia_input,6,'FaceColor','k','FaceAlpha',0.3,'Orientation','horizontal');axis square;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);xlim([0 70]);
 %% Display Ex and IN and difference sorted by pial depth
display_sortfr([pia_input L5fr],1,'L5')
%% Plot fraction vs pia
plot_fraction(L23fr,L4fr,L5fr,pia_input);
%% Plot average maps all
 plot_avg_maps(str,1:147,ex_map,in_map,pia_input,10,0,[]);
 
 %% Discretize maps based on pia
[gr gg]=discretize(pia_input,3);
plot_avg_maps(str,find(gr==1),ex_map,in_map,pia_input,10,0,[]);
plot_avg_maps(str,find(gr==2),ex_map,in_map,pia_input,10,0,[]);
plot_avg_maps(str,find(gr==3),ex_map,in_map,pia_input,10,0,[]);
%% in vivo and maps

%% Creating correlation matrix
com=[];com=[L23fr(:,1)  L23fr(:,2) L4fr(:,1)  L4fr(:,2) L5fr(:,1) L5fr(:,2) diffL23fr diffL4fr diffL5fr abs(od_out_iviv(:,3)) pia_input]; 
G=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);

%% PCI vs L4EX-IN
par1=com(:,14)
kk1=find(isnan(par1))
par2=com(:,8)
par3=[str(:).resp]
kk2=find(isnan(par2))
par1(kk1)=[];
par2(kk1)=[];
par3(kk1)=[];
corr_plot(par1,par2,[],{'PCI','L4fr EX-IN',''})
%% %% Circular correlation for ORI
a=find(od_out_iviv(:,1)>0.25); 
par_c=[];
 par_c=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1)  L5fr(a,2) diffL23fr(a) diffL4fr(a) diffL5fr(a) ]
for i=1:9
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(a,4)), par_c(:,i))
end
% Circular correlation for DRI
a=[];
a=find(od_out_iviv(:,2)>0.25); 
par_c=[];
par_c=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1)  L5fr(a,2) diffL23fr(a) diffL4fr(a) diffL5fr(a) ]
for i=1:9
    [rho2(i) pval2(i)] = circ_corrcl(deg2rad(od_out_iviv(a,5)), par_c(:,i))
end
%% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% 
a=[];
a=find(od_out_iviv(:,1)>0); 
figure;scatter(diffL5fr(a,1),abs(od_out_iviv(a,3)),20,od_out_iviv(a,2),'filled')
%% 
a=[];
a=find(od_out_iviv(:,1)>0); 
figure;scatter(L23fr(a,1),abs(od_out_iviv(a,3)),20,pia_input(a),'filled')
%% 
eye_sp=[contra_id ipsi_id];

par=diffL5fr
[statsout]=dual_barplot(par,eye_sp,bino_id,0);xticks([1:1:2]);hold on;xticklabels({'mono' ,'bino'});ylabel('EX-IN L5');set(gca,'FontSize',10);xtickangle(45);
%% 
eye_sp=[contra_id ipsi_id];

par=diffL23fr
[statsout]=dual_barplot(par,eye_sp,bino_id,0);xticks([1:1:2]);hold on;xticklabels({'mono' ,'bino'});ylabel('EX-IN L5');set(gca,'FontSize',10);xtickangle(45);

%% 
par=[]
par=morph_parameters(:,19)
%par=ang_inL5(:,8)
%par=span(:,3)*69
par=pia_input
[statsout]=dual_barplot(par,eye_sp,bino_id,0);xticks([1:1:2]);hold on;xticklabels({'mono' ,'bino'});ylabel('EX-IN L5');set(gca,'FontSize',10);xtickangle(45);

%% 
corr_plot(diffL23fr,diffL5fr,[],{'PCI','L4fr EX-IN',''})


%% 
figure
s1=eye_sp
frac_h=[];
frac_h=frv(:,:);
frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
%EX and IN Horizontal
;hold on;
% for i=1:size(frac_h(s1,:),1)
% exp=plot(frac_h(s1(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s1(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
mexp=errorbar(nanmean(frac_h(s1,1:16)),nanstd(frac_h(s1,1:16))/sqrt(size(frac_h(s1,:),1)),'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;

hold on;
mexp=errorbar(nanmean(frac_h(s1,17:end))*-1,nanstd(frac_h(s1,17:end))/sqrt(size(frac_h(s1,:),1))*-1,'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
hold on;
s2=bino_id
frac_h=[];
frac_h=frv(:,:);
frac_diffh=frac_h(s2,1:16)-frac_h(s2,17:end);
% for i=1:size(frac_h(s2,:),1)
% exp=plot(frac_h(s2(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s2(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
hold on;
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s2,:),1)),'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);

hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'k');
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














