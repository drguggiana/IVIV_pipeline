%% START
%Load structure with 147 cells used for the paper
str      = 'D:\Postdoc_Margrie\Projects\L23\structure';
folder_list = uipickfiles('FilterSpec',str);
load(char(folder_list));
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
%Get the morphology traces
for i=1:length(str)
    if ~isempty(str(i).morph)==1 
    morph_traces{:,i}=str(i).morphtraces;
    else
        morph_traces{:,i}=NaN;
    end
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
 %total sum normlaized for ex and inh 
for i=1:length(str)
%whole map
tmp1=ex_map(3:end,:,i);
tmp2=in_map(:,:,i);
% tmp1=ex_map(:,:,i);
% tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
%tot_input(i,:)=[sum(tmp1(:))/length(tmp1);sum(tmp2(:))/length(tmp2)];
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map(3:5,:,i);
tmp2=in_map(3:5,:,i);
tot_inputL23(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map(6:7,:,i);
tmp2=in_map(6:7,:,i);
tot_inputL4(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L5
tmp1=ex_map(8:11,:,i);
tmp2=in_map(8:11,:,i);
tot_inputL5(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
end
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
r2=1-(sum(sqrt((y-yp_s').^2))/sum(sqrt((y-ym).^2)));
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
%%  Histogram of pia distribution with color coded in vivo morpho by itself
close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 300, 300]);
hold on;histogram(pia_input(iviv_cells),'FaceColor','k','EdgeColor','k','Orientation','horizontal');
hold on;histogram(pia_input(find([str(:).iviv]==1 & morph_cells==1)),'FaceColor','w','EdgeColor','k','Orientation','horizontal');
legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
    [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
ylim([100 350])
legend boxoff ;set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);
;ylabel('Pial depth (�m)');xlabel('Cell count');box off;
%% FIGURE 2
%% Show vertical per layer for iviv cells
 [stats_g] = display_inputs_part2([frv(iviv_cells,:)],[frh(iviv_cells,:)],frv(iviv_cells,1:16)-frv(iviv_cells,17:end),frh(iviv_cells,1:16)-frh(iviv_cells,17:end),[]);
%%  Plot average iviv cells maps
plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
%% Plot example maps
cnr=1;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);

%%  Extract horizontal fraction per layer
[L1h L23h L4h L5h] = h_fraclayer(ex_map, in_map);
%% Plot horizontal fraction per layer
plot_horizontal_fraction(iviv_cells,L23h,L4h,L5h);

%% Plot horizontal centroid per layer with ref line and exclude outliers
g=find(ang_inL4(iviv_cells,3)*69-ang_inL4(iviv_cells,1)*69<135)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
plot(ang_exL23(iviv_cells,3)*69-ang_exL23(iviv_cells,1)*69,ang_inL23(iviv_cells,3)*69-ang_inL23(iviv_cells,1)*69,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-160 160];xticks([-160:80:160]);yticks([-160:80:160]);
ref.YData=[-160 160];
ref.Color='k';box off;xlim([-160 160]);ylim([-160 160]);hold on;title('L23','FontWeight','normal');%xticks([-100:50:100]);yticks([-100:50:100]);
subplot(1,3,2);
plot(ang_exL4(iviv_cells,3)*69-ang_exL4(iviv_cells,1)*69,ang_inL4(iviv_cells,3)*69-ang_inL4(iviv_cells,1)*69,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-160 160];
ref.YData=[-160 160];ylabel('IN C_{x} (�m)','Color','b')
ref.Color='k';box off;xlim([-160 160]);ylim([-160 160]);hold on;title('L4','FontWeight','normal');xticks([-160:80:160]);yticks([-160:80:160]);
xlabel('EX C_{x} (�m)','Color','r')
subplot(1,3,3);
%L5, remove outlier
tempL5ex=ang_exL5(iviv_cells,3)*69-ang_exL5(iviv_cells,1)*69;
tempL5in=ang_inL5(iviv_cells,3)*69-ang_inL5(iviv_cells,1)*69;
tempL5ex(find(tempL5ex<-400))=[];
tempL5in(find(tempL5in<-400))=[];
plot(tempL5ex,tempL5in,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;;hold on;title('L5','FontWeight','normal');%xticks([-300:150:300]);yticks([-300:150:300]);
xlim([-160 160]);ylim([-160 160]);xticks([-160:80:160]);yticks([-160:80:160]);
 ref.XData=[-160 160];
 ref.YData=[-160 160];
 %% Orientation tuning
 %% Plot correlations between orientation selectivity/tuning width and translaiminr inputs
 %Fraction 
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(L23fr(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(L23fr(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);ylabel('gOSI');
hold on;text(0.4,1,'L2/3');set(gca,'FontSize',10);
subplot(1,3,2);
scatter(L4fr(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(L4fr(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);xlabel('Fraction of total input');
hold on;text(0.3,1,'L4');set(gca,'FontSize',10);xlim([-0.02 0.6])
subplot(1,3,3);
scatter(L5fr(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(L5fr(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);xlim([-0.02 0.6])
hold on;text(0.3,1,'L5');set(gca,'FontSize',10);legend('IN','EX');legend box off
%% 
%Lateral extent
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(span(:,4),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(span(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);ylabel('gOSI');
subplot(1,3,2);
scatter(span(:,5),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(span(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);xlabel('Horizontal extent');
subplot(1,3,3);
scatter(span(:,6),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(span(:,3),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);
%Total count of normalied input
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
scatter(tot_inputL23(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL23(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);ylabel('gOSI');
subplot(1,3,2);
scatter(tot_inputL4(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL4(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);xlabel('Horizontal extent');
subplot(1,3,3);
scatter(tot_inputL5(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL5(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);
%% Group selective and unselective cells
g1=[];g2=[];
t1=find(od_out_iviv(:,1)<=0.25)
g2=find(od_out_iviv(:,1)>0.25)
g1=[unres_id t1'];
par=L5fr(:,2);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);hold on;xticklabels({'UT' ,'TU'});ylabel('C_{x} (�m)');set(gca,'FontSize',10);xtickangle(45);
%% Plot vertical profile selective vs unselective profiles
frac_v=frv;
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [300,100, 250, 300]);set(gcf,'color','w');
mexp=errorbar(nanmean(frac_v(g1,1:16)),1:16,nanstd(frac_v(g1,1:16))/sqrt(size(frac_v(g1,1:16),1)),'horizontal','r');set(gca,'Ydir','reverse');
mexp.CapSize=3;
% Set transparency level (0:1)
alpha = 0.3;   
% Set transparency (undocumented)
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha])
 hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
 hold on;line([-1 1], [9 9],'Color','k','LineStyle','--');
mexp2=errorbar(nanmean(frac_v(g1,17:end))*-1,1:16,nanstd(frac_v(g1,17:end))/sqrt(size(frac_v(g1,17:end),1)),'horizontal','b');set(gca,'Ydir','reverse');
% Set transparency level (0:1)
alpha = 0.3;   
% Set transparency (undocumented)
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha])
mexp2.CapSize=3;xlim([-0.4 0.4]);xticks([-0.4:0.2:0.4]);hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});xlabel('Fraction of total input');set(gca,'FontSize',10);ylabel('Depth (�m)');set(gca,'FontSize',10)
hold on;
mexp=errorbar(nanmean(frac_v(g2,1:16)),1:16,nanstd(frac_v(g2,1:16))/sqrt(size(frac_v(g2,1:16),1)),'horizontal','r','Linestyle', '--');set(gca,'Ydir','reverse');
 hold on;line([-1 1], [3 3],'Color','k','LineStyle','--');hold on;line([-1 1], [6 6],'Color','k','LineStyle','--');
 hold on;line([-1 1], [9 9],'Color','k','LineStyle','--');
 alpha = 0.8;   
% Set transparency (undocumented)
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha])
 mexp.CapSize=3;
mexp2=errorbar(nanmean(frac_v(g2,17:end))*-1,1:16,nanstd(frac_v(g2,17:end))/sqrt(size(frac_v(g2,17:end),1)),'horizontal','b','Linestyle', '--');set(gca,'Ydir','reverse');
mexp2.CapSize=3;
alpha = 0.8;   
% Set transparency (undocumented)
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha]);
box off;

%% 



%% Direction tuning asymmetry

%% Rolling average
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (�m)','Color','k');hold on;xlabel('Direction preference (deg)')
hold on;set(gca,'FontSize',10)
%% Bar plot Cx for IN and EX in L23
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
 binodir_delta=binodir(:,1)-binodir(:,2);
% 
idx_bi=find(abs(od_out_iviv(a,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
a(idx_bi(id_ov))=[];
% sector=60;
% midpoint=130;
 sector=60;
 midpoint=130;
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
g1=[];
g2=[];
g3=[];
g4=[];
g1=find(od_out_iviv(a,5)>s1a & od_out_iviv(a,5)<s1b) ;
g2=find(od_out_iviv(a,5)>s2a & od_out_iviv(a,5)<s2b);
g3=find(od_out_iviv(a,5)>s3a & od_out_iviv(a,5)<s3b);
g4=find(od_out_iviv(a,5)>s4a & od_out_iviv(a,5)<s4b);
par=[];
s1=[g1' g2'];
s2=[g3' g4'];
%par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (�m)');set(gca,'FontSize',10);xtickangle(45);ylim([-70 70]);yticks(-70:35:70)
%% Plot fraction for aligned and non aligned 
statsout=plot_horizontal_fraction_group(str,od_out_iviv,L23h,L4h,L5h)




%% Morphology horizontal sholl analysis from VS
%align all morhologies correctly using m_flip again (morph_flip is just a
%147 matlab table with 1, 0 and NaN)
cd('D:\Postdoc_Margrie\Projects\L23\structure\structure_with_morphflip info');
load('morph_flip.mat');
%% load cyl_coord
cd('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Manuscript\split_paper\short_paper\Volker_xy_coord_branchpoint');
load('cyl_coord.mat');
%% orient all morpholgies medial to lateral
for i=1:length(str)
    if morph_flip(i)==1
mx_basal{:,i}= [cyl_coord.basal{1, i}(:,1:2)*-1 cyl_coord.basal{1, i}(:,3:end)];
mx_apical{:,i}= [cyl_coord.apical{1, i}(:,1:2)*-1 cyl_coord.apical{1, i}(:,3:end)];
    else
        mx_basal{:,i}= cyl_coord.basal{1, i};
        mx_apical{:,i}= cyl_coord.apical{1, i};
    end
end
%% Plot examples with to check flip of morphology
figure;
cr=1;scatter(cyl_coord.basal{1, cr}(:,1),cyl_coord.basal{1, cr}(:,3));
hold on;scatter(cyl_coord.apical{1, cr}(:,1),cyl_coord.apical{1, cr}(:,3));
%% 
cr=3;
figure;scatter(mx_basal{1, cr}(:,1),mx_basal{1, cr}(:,3));
hold on;scatter(mx_apical{1, cr}(:,1),mx_apical{1, cr}(:,3));
set(gcf,'color','w');box off;
axis off;set(gca, 'YDir','reverse');
data=mx_basal{1, cr}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
  figure;plot(output.resampleX,output.resampleScholl);
  set(gcf,'color','w');box off;ylabel('Horizontal sholl: Branch points');
  data=mx_apical{1, cr}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
 hold on; plot(output.resampleX,output.resampleScholl);
%% 
%% Read out desired cells for direction part
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
 binodir_delta=binodir(:,1)-binodir(:,2);
idx_bi=find(abs(od_out_iviv(a,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
a(idx_bi(id_ov))=[];
 sector=60;
 midpoint=130;
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
g1=[];
g2=[];
g3=[];
g4=[];
g1=find(od_out_iviv(a,5)>s1a & od_out_iviv(a,5)<s1b) ;
g2=find(od_out_iviv(a,5)>s2a & od_out_iviv(a,5)<s2b);
g3=find(od_out_iviv(a,5)>s3a & od_out_iviv(a,5)<s3b);
g4=find(od_out_iviv(a,5)>s4a & od_out_iviv(a,5)<s4b);
par=[];
s1=[g1' g2'];
s2=[g3' g4'];
%% Do some further extraction and direction aligning 
for i=1:length(a)
te{:,i}= mx_apical{1, a(i)};
te2{:,i}= mx_basal{1, a(i)};
end
%for groups s1 and s2
for i=1:length(s1)
s1_te{:,i}=te{1,s1(i)};
s1_te2{:,i}=te2{1,s1(i)};
end
%flip acccording to direction 315
s1_te{:,1}=s1_te{:,1}*-1;
s1_te{:,2}=s1_te{:,2}*-1;

for i=1:length(s2)
s2_te{:,i}=te{1,s2(i)};
s2_te2{:,i}=te2{1,s2(i)};
end
s1_coord.apical=s1_te;
s1_coord.basal=s1_te2;
s2_coord.apical=s2_te;
s2_coord.basal=s2_te2;
%% Plot aligned and non aligned groups
s1_out=shell_scholl_analysis_horizontal_vs1(s1_coord,5,300);
s2_out=shell_scholl_analysis_horizontal_vs1(s2_coord,5,300);
close all;
%% Direction biases
shadedErrorBar(s1_out.resampleX,s1_out.meanScholl,s1_out.stdScholl,'lineProps','m');
hold on;
shadedErrorBar(s2_out.resampleX,s2_out.meanScholl,s2_out.stdScholl,'lineProps','g');
legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Basal tree')
%title('Apical tree')
legend('Aligned','NonAligned');ylabel('Horizontal sholl: Branch points');

%
%% sum medial and lateral
sampling=5;
extend=300;
for i=1:length(s1)

    data=[];
    data=s1_coord.basal{:,i}(:,1:2);
s1_outdir=scholl_analysis_horizontal_vs1(data,sampling,extend)
temp=[sum(s1_outdir.resampleScholl(find(s1_outdir.resampleX<0 & s1_outdir.resampleX>-50))) sum(s1_outdir.resampleScholl(find(s1_outdir.resampleX<50 & s1_outdir.resampleX>0)))];
ml_basal(:,i)=temp;
end
stats=paired_plot(ml_basal',0,{'r','k'});


%% Read out orientation for morphology
a=[];
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0); 
sector=45;
midpoint=130;
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint)-sector/2];s4b=[(midpoint)+sector/2];
g1=[];
g2=[];
g3=[];
g4=[];
s1=[];
s2=[];
g3=find(od_out_iviv(a,4)>s3a & od_out_iviv(a,4)<s3b);
g4=find(od_out_iviv(a,4)>s4a & od_out_iviv(a,4)<s4b);
par=[];
s1=[g3'];
s2=[g4'];
%% Orientation 
for i=1:length(a)
te_or{:,i}= mx_apical{1, a(i)};
te2_or{:,i}= mx_basal{1, a(i)};
end
%for groups s1 and s2
for i=1:length(s1)
s1_te_or{:,i}=te_or{1,s1(i)};
s1_te2_or{:,i}=te2_or{1,s1(i)};
end

for i=1:length(s2)
s2_te_or{:,i}=te_or{1,s2(i)};
s2_te2_or{:,i}=te2_or{1,s2(i)};
end
s1_coord_or.apical=s1_te_or;
s1_coord_or.basal=s1_te2_or;
s2_coord_or.apical=s2_te_or;
s2_coord_or.basal=s2_te2_or;
%% Plot aligned and non aligned groups
s1_out_or=shell_scholl_analysis_horizontal_vs1(s1_coord_or,5,300)
s2_out_or=shell_scholl_analysis_horizontal_vs1(s2_coord_or,5,300);
close all;
%% Plot averages with shaded errorbar
shadedErrorBar(s1_out_or.resampleX,s1_out_or.meanScholl,s1_out_or.stdScholl,'lineProps','m');
hold on;
shadedErrorBar(s2_out_or.resampleX,s2_out_or.meanScholl,s2_out_or.stdScholl,'lineProps','g');
legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Basal tree')
title('Apical tree')
legend('Aligned','NonAligned');ylabel('Horizontal sholl: Branch points');
%% Test pial depth 
s1=[];
s2=[];
s1=find(pia_input<170);
s2=find(pia_input>170);

for i=1:length(s1)
s1_te_pia{:,i}= mx_apical{1, s1(i)};
s1_te2_pia{:,i}= mx_basal{1, s1(i)};
end

for i=1:length(s2)
s2_te_pia{:,i}= mx_apical{1, s2(i)};
s2_te2_pia{:,i}= mx_basal{1, s2(i)};
end
s1_coord_pia.apical=s1_te_pia;
s1_coord_pia.basal=s1_te2_pia;
s2_coord_pia.apical=s2_te_pia;
s2_coord_pia.basal=s2_te2_pia;
%% Plot upper vs lower
s1_out_pia=shell_scholl_analysis_horizontal_vs1(s1_coord_pia,5,300)
s2_out_pia=shell_scholl_analysis_horizontal_vs1(s2_coord_pia,5,300);
close all;
%% Plot averages with shaded errorbar
shadedErrorBar(s1_out_pia.resampleX,s1_out_pia.meanScholl,s1_out_pia.stdScholl,'lineProps','m');
hold on;
shadedErrorBar(s2_out_pia.resampleX,s2_out_pia.meanScholl,s2_out_pia.stdScholl,'lineProps','g');
legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Basal tree')
title('Apical tree')
legend('Upper','Lower');ylabel('Horizontal sholl: Branch points');





%% Supplements
close all;
for i=1:length(iviv_cells)
figure;
iviv_plotter(str,iviv_cells(i),i)
end
%% Save individual morpho, in vivo, maps panels
fn='D:\Postdoc_Margrie\Projects\L23\Paper1\Supp1'
savepdf_SW(fn,0);
%% Supp Fig 4: Cell to cell fraction and horzontal span 
%% Plot difference between vertical excitation and inhbition per cell
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
subplot(3,1,1);
histogram(diffL23fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL23fr),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L2/3');ylabel('Cell Counts');set(gca,'FontSize',10);text(-0.5,28,'IN','Color','b');
set(gca,'FontSize',10);text(0.5,28,'EX','Color','r');set(gca,'FontSize',10)
subplot(3,1,2);
histogram(diffL4fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL4fr),30,'kv','filled');ylabel('Cell Counts')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L4');set(gca,'FontSize',10);set(gca,'FontSize',10);
subplot(3,1,3);
histogram(diffL5fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL5fr),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L5');set(gca,'FontSize',10);ylabel('Cell Counts');xlabel('EX-IN vertical fraction')
%stats against 0
[p,h] = signrank(diffL5fr);
%% %% Plot difference between horizintal excitation and inhbition per cell SPAN
diffL23fr_span=span(:,1)-span(:,4)
diffL4fr_span=span(:,2)-span(:,5)
diffL5fr_span=span(:,3)-span(:,6)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 600]);
subplot(3,1,1);
histogram(diffL23fr_span(iviv_cells)*69,'BinWidth',69,'FaceColor','w');box off;ylim([0 30]);
xlim([-600 600]);xticks([-600:300:600]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL23fr_span*69),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-500,24,'L2/3');ylabel('Cell Counts');set(gca,'FontSize',10)
subplot(3,1,2);
histogram(diffL4fr_span(iviv_cells)*69,'BinWidth',69,'FaceColor','w');box off;ylim([0 30]);
xlim([-600 600]);xticks([-600:300:600]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL4fr_span*69),30,'kv','filled');ylabel('Cell Counts')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-500,24,'L4');set(gca,'FontSize',10);;set(gca,'FontSize',10);
text(-300,28,'IN','Color','b');set(gca,'FontSize',10);text(300,28,'EX','Color','r');set(gca,'FontSize',10)
subplot(3,1,3);
histogram(diffL5fr_span(iviv_cells)*69,'BinWidth',69,'FaceColor','w');box off;ylim([0 30]);
xlim([-600 600]);xticks([-600:300:600]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL5fr_span*69),30,'kv','filled');ylabel('Cell Counts')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-500,24,'L5');set(gca,'FontSize',10);
xlabel('EX-IN horizontal extent (�m)');
%stats against 0
[p,h] = signrank(diffL5fr_span);