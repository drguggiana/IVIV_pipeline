
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
;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
 %% Plot average iviv cells maps
 plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
 %% Plot example maps
  plot_avg_maps(str,116,ex_map,in_map,pia_input,1,0,[]);
%% Show vertical per layer for iviv cells
 [stats_g] = display_inputs_part2([frv(iviv_cells,:)],[frh(iviv_cells,:)],frv(iviv_cells,1:16)-frv(iviv_cells,17:end),frh(iviv_cells,1:16)-frh(iviv_cells,17:end),[]);
  %% Show horizontal fraction per layer
[L1h L23h L4h L5h] = h_fraclayer(ex_map, in_map);
%% Plot horizontal fraction per layer
frac_h=L23h(iviv_cells,:);
frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position')
hold on;text(2.5,0.6,'L2/3','FontSize',12);
frac_h=L4h(iviv_cells,:);hold on;

frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
subplot(1,3,2);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');%ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');%ylabel('Fraction of total input');
box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position');
hold on;text(2.5,0.6,'L4','FontSize',12);

frac_h=L5h(iviv_cells,:);

frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
subplot(1,3,3);hold on;
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');%ylabel('Fraction of total input');
box off;
exp.Color(4) = 0.05;
exp=plot(frac_h(i,17:end)'*-1,'-b');%ylabel('Fraction of total input');box off;
exp.Color(4) = 0.05;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-0.6 0.6]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
set(gca,'FontSize',10);
hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position')
hold on;text(2.5,0.6,'L5','FontSize',12);

%% Plot centroid Cx fdr EX IN
%% L23 and L4 centroid EX and IN
corr_plot((ang_exL23(iviv_cells,3)-ang_exL23(iviv_cells,1))*69,(ang_inL23(iviv_cells,3)-ang_inL23(iviv_cells,1))*69,[],{'EX','IN',''});ylabel('IN C_{x} (µm)','Color','b');xlabel('EX C_{x} (µm)','Color','r');
temp1=ang_exL4(iviv_cells,3)-ang_exL4(iviv_cells,1);
temp2=ang_inL4(iviv_cells,3)-ang_inL4(iviv_cells,1);
temp1(find(isnan(temp2)))=[];
temp2(find(isnan(temp2)))=[];
corr_plot(temp1*69,temp2*69,[],{'EX','IN',''});ylabel('IN C_{x} (µm)','Color','b');xlabel('EX C_{x} (µm)','Color','r');
%% L5 centroid
temp1=[];temp2=[];
temp1=ang_exL5(iviv_cells,3)-ang_exL5(iviv_cells,1);
temp2=ang_inL5(iviv_cells,3)-ang_inL5(iviv_cells,1);
aa=[find(isnan(temp2)); find(isnan(temp1))]
temp1(aa)=[];
temp2(aa)=[];
corr_plot(temp1*69,temp2*69,[],{'EX','IN',''});ylabel('IN','Color','b');xlabel('EX','Color','r');%xlim([-200 200]);ylim([-200 200]);xlabel('EX C_{x} (µm)','Color','r');
%% Plot horizontal centroid per layer with ref line 
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
ref.YData=[-160 160];ylabel('IN C_{x} (µm)','Color','b')
ref.Color='k';box off;xlim([-160 160]);ylim([-160 160]);hold on;title('L4','FontWeight','normal');xticks([-160:80:160]);yticks([-160:80:160]);
xlabel('EX C_{x} (µm)','Color','r')
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
%% Rolling average
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm)','Color','k');hold on;xlabel('Direction preference (deg)')
hold on;set(gca,'FontSize',10)
%% Circular correlation 
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_inL23(a,3)-ang_inL23(a,1))*69))
%% ipsi contra
tmp = cat(3,[str(:).contra],[str(:).ipsi]); 
eye_sp = nansum(tmp,3);
tmp2 = cat(3,eye_sp,[str(:).bino]*2);
eye_spec=nansum(tmp2,3);
%% bino cells
bino_dir=reshape([str(find(eye_spec==2)).Dir],[20,2])
stats=paired_plot(bino_dir,1,{'b','r'})
%% ODI of cells used for comparison of driections
corr_plot(od_out_iviv(a,4)-45,abs((ang_inL23(a,3)-ang_inL23(a,1))*69),od_out_iviv(a,3),{'Orientation along slice','IN',''});
%% Absolute Cx relative to axis of displacement
close all;
%a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
 binodir_delta=binodir(:,1)-binodir(:,2);
% 
idx_bi=find(abs(od_out_iviv(a,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
a(idx_bi(id_ov))=[];
fig3= figure;set(fig3, 'Name', 'Input distribution');set(fig3, 'Position', [200, 600, 200, 200]);set(gcf,'color','w');
scatter(od_out_iviv(a,4)-45,abs((ang_inL23(a,3)-ang_inL23(a,1))*69),12,'bo','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
%hold on;scatter(od_out_iviv(a,4)-45,(ang_inL23(a,3)-ang_inL23(a,1))*69,12,'ro','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_inL23(a,3)-ang_inL23(a,1))*69));
text(-20,90,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','b')
hold on
scatter(od_out_iviv(a,4)-45,abs((ang_exL23(a,3)-ang_exL23(a,1))*69),12,'ro','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_exL23(a,3)-ang_exL23(a,1))*69));
text(-20,100,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','r')
set(gca,'FontSize',10);
ylim([-10 100])
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
par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% 
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
b=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)<0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
 binodir_delta=binodir(:,1)-binodir(:,2);
% 
idx_bi=find(abs(od_out_iviv(a,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
a(idx_bi(id_ov))=[];

binodir_b=reshape([str(b).Dir],[2,length(a)])';
 binodir_delta_b=binodir_b(:,1)-binodir_b(:,2);
% 
idx_bi_b=find(abs(od_out_iviv(b,3))<0.1)
id_ov_b=find(abs(binodir_delta_b(find(abs(od_out_iviv(b,3))<0.1)))>90)
b(idx_bi_b(id_ov_b))=[]
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
g3=find(od_out_iviv(b,5)>s1a & od_out_iviv(b,5)<s1b) ;
g4=find(od_out_iviv(b,5)>s2a & od_out_iviv(b,5)<s2b);
par=[];
s1=[g1' g2'];
s2=[g3' g4'];
par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g1)=par(g1)*-1
par(g3)=par(g3)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% 
par2=[];
%a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0); 
par2=od_out_iviv(a,2);
corr_plot(abs(par(s1)),par2(s1),[],{'Cx offset','DSI','Pial depth'});

%% Bar plot 
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0); 
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
par=eye_spec(a)
[statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% Correlation with DSI
a=[];
par2=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0); 
par2=od_out_iviv(a,2);
corr_plot(abs(par(s1)),par2(s1),[],{'Cx offset','DSI','Pial depth'});
%% EX and IN
par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
par1(g1)=par1(g1)*-1;
par2=(ang_exL23(a,3)-ang_exL23(a,1))*69;
par2(g1)=par2(g1)*-1;
p1=paired_plot([par2(s1) par1(s1)],1,{'r','b'});xticklabels({'EX','IN'});ylabel('C_{x} (µm)');set(gca,'FontSize',10)

%% Correlation comparison across layers
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
g1=[];g2=[];
g1=s1;
g2=s2;
par1=[];par2=[];
par1=ang_exL23(a,3)*69-ang_exL23(a,1)*69
par2=ang_exL4(a,3)*69-ang_exL4(a,1)*69
[R1,P1,RLO1,RUP1]= corrcoef(par1(g1), par2(g1),'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_inL23(a,3)*69-ang_inL23(a,1)*69
par2=ang_inL4(a,3)*69-ang_inL4(a,1)*69
[R4,P4,RLO4,RUP4]= corrcoef(par1(g1), par2(g1),'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_exL23(a,3)*69-ang_exL23(a,1)*69
par2=ang_exL4(a,3)*69-ang_exL4(a,1)*69
[R7,P7,RLO7,RUP7]= corrcoef(par1(g2), par2(g2),'Rows','pairwise', 'alpha', 0.05);
par1=[];par2=[];
par1=ang_inL23(a,3)*69-ang_inL23(a,1)*69
par2=ang_inL4(a,3)*69-ang_inL4(a,1)*69
[R10,P10,RLO10,RUP10]= corrcoef(par1(g2), par2(g2),'Rows','pairwise', 'alpha', 0.05);
%Figure
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 200, 300]);
scatter([R1(2); R4(2)],[1 3],'m','filled')
hold on;scatter([R7(2); R10(2)],[2 4],'c','filled')
ylim([0 5]);yticks([1:1:5]);xlabel('Correlation');
hold on;plot([RLO1(2) RUP1(2)],[1 1],'-m');
hold on;plot([RLO4(2) RUP4(2)],[3 3],'-m');
hold on;plot([RLO7(2) RUP7(2)],[2 2],'-c');
hold on;plot([RLO10(2) RUP10(2)],[4 4],'-c');
yticklabels({'EX','EX','IN','IN'});title('C_{x} L2/3 to C_{x} L4');
legend('AL','NAL');legend boxoff;
set(gca,'FontSize',11);

%% Plot horizontal fraction per layer
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
%START
frac_h=[];
frac_h=L23h(a,:);
flip_hl=[flip(frac_h(g1,1:16),2); frac_h(g2,1:16)];
flip_hl_in=[flip(frac_h(g1,17:end),2); frac_h(g2,17:end)]
frac_diffh=[];
% frac_diffh=frac_h(flip_hl,1:16)-frac_h(s1,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 600, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
% for i=1:size(frac_h(s1,:),1)
% exp=plot(frac_h(s1(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s1(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
%mexp=errorbar(nanmean(frac_h(s1,1:16)),nanstd(frac_h(s1,1:16))/sqrt(size(frac_h(s1,:),1)),'Color',[0.7 0 0.4]);
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
% hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
% mexp.CapSize=3;
ylabel('Fraction of total input');
frac_h=[];
frac_h=L23h(a,:);
frac_diffh=[];
frac_diffh=frac_h(s2,1:16)-frac_h(s2,17:end);
%EX and IN Horizontal
subplot(1,3,1);hold on;
% for i=1:size(frac_h(s2,:),1)
% exp=plot(frac_h(s2(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s2(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
hold on;
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s1,:),1)),'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
%hold on;text(2.5,0.4,'L2/3','FontSize',12);
% hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh(:,:),1)),'b');
% mexp.CapSize=3;


%L4
frac_h=[];
frac_h=L4h(a,:);
flip_hl=[flip(frac_h(g1,1:16),2); frac_h(g2,1:16)];
flip_hl_in=[flip(frac_h(g1,17:end),2); frac_h(g2,17:end)]
%frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
%EX and IN Horizontal
subplot(1,3,2);hold on;
% for i=1:size(frac_h(s1,:),1)
% exp=plot(frac_h(s1(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s1(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
frac_h=[];
frac_h=L4h(a,:);
frac_diffh=frac_h(s2,1:16)-frac_h(s2,17:end);
% for i=1:size(frac_h(s2,:),1)
% exp=plot(frac_h(s2(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s2(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
hold on;
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s1,:),1)),'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'Color',[0 0.5 0.5]);
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
%hold on;text(2.5,0.55,'L4','FontSize',12);
xlabel('Horizontal map position')

%L5
frac_h=[];
frac_h=L5h(a,:);
flip_hl=[flip(frac_h(g1,1:16),2); frac_h(g2,1:16)];
flip_hl_in=[flip(frac_h(g1,17:end),2); frac_h(g2,17:end)]
%frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
%EX and IN Horizontal
subplot(1,3,3);hold on;
% for i=1:size(frac_h(s1,:),1)
% exp=plot(frac_h(s1(i),1:16)','-r');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% exp=plot(frac_h(s1(i),17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
% exp.Color(4) = 0.05;
% end
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color',[0.7 0 0.4]);
hold on;mexp.CapSize=3;
frac_h=[];
frac_h=L5h(a,:);
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
%hold on;text(2.5,0.55,'L5','FontSize',12);
%% Statistics
par=[];
par=L23h(a,10);

if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par(s1),par(s2));
    statsout=p;
elseif t==1
    [p k]=ranksum(par(s1),par(s2),'tail','right');
    statsout=p;
else t==2
    [p k]=ranksum(par(s1),par(s2),'tail','left');
    statsout=p;
end


%% 

%% Cecking all morpho parameters if needed
miviv=find([str(:).iviv]==1 & morph_cells==1);
smo=morph_cells(a)
par=[];
par=morph_parameters(a,1)
[statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);



%% Looking at span
parameter_vector = span(:,1)*69;
parameter_vector2 = span(:,4)*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1);
ylabel('Horizontal extent')
%% Span barplot
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
par=span(a,1)*69;
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})
%% 


%% Check visual non visual cells
sector=60;
midpoint=130;
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25 & od_out_iviv(:,5)>s1a & od_out_iviv(:,5)<s1b); 
b=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25 & od_out_iviv(:,5)>s2a & od_out_iviv(:,5)<s2b)
s1=[a' b'];
s2=unres_id;
par=(ang_inL23(:,3)-ang_inL23(:,1))*69;
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})
%% 
sector=60;
midpoint=130;
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25 & od_out_iviv(:,5)>s3a & od_out_iviv(:,5)<s3b); 
b=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25 & od_out_iviv(:,5)>s4a & od_out_iviv(:,5)<s4b)
s1=[a' b'];
s2=unres_id;
par=abs((ang_inL23(:,3)-ang_inL23(:,1))*69);
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})

%% 
a=[];
b=[];
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0.25); 
b=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25 & od_out_iviv(:,5)>s2a & od_out_iviv(:,5)<s2b)
s1=[a'];
s2=unres_id;
par=L1fr(:,2);
[statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})

%% Using asymmerty as calculated by Volker S
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0 & abs(od_out_iviv(:,3))>0); 
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

frac_h=L23h(:,:);

for i=1:length(frac_h)
    frac_asy(i)=asymmetry_quantification(frac_h(i,1:16))
    frac_asy_in(i)=asymmetry_quantification(frac_h(i,17:end))
end

par=[frac_asy(a).asymmetricComponent]-[frac_asy_in(a).asymmetricComponent]


[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})
%% 
%% EX and IN
% par1=(ang_exL23(a,3)-ang_exL23(a,1))*69;
% par2=(ang_inL5(a,3)-ang_inL5(a,1))*69;
% par2(g1)=par2(g1)*-1;
% par1(g1)=par1(g1)*-1;
par1=span(a,4)
par2=span(a,1)
p1=paired_plot([par2(s1) par1(s1)],1,{'r','b'});xticklabels({'EX','IN'});ylabel('C_{x} (µm)');set(gca,'FontSize',10)

%% Correlation with DSI
par=abs(span(a,1)-span(a,4))
par2=[];
%a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0); 
par2=od_out_iviv(a,2);
corr_plot(abs(par(s1))',par2(s1),[],{'Cx offset','DSI','Pial depth'});
%% Morphology
close all;
for i=1:length(str)
    if ~isempty(str(i).morph)==1 
    zz{:,i}=str(i).morphtraces;
    else
        zz{:,i}=NaN;
    end
end
%% Basal
for i=1:147
   
X_coord=[];
Y_coord=[];
try
[B X_coord Y_coord]=B_tree_sw (zz{:,i}{2}, '-s');
x_all_basal{:,i}=X_coord;
y_all_basal{:,i}=Y_coord;
catch
  x_all_basal{:,i}=NaN;
y_all_basal{:,i}=NaN;  
end
end
%% Apical
for i=1:147
   
X_coord=[];
Y_coord=[];
try
[B X_coord Y_coord]=B_tree_sw (zz{:,i}{1}, '-s');
x_all_apical{:,i}=X_coord;
y_all_apical{:,i}=Y_coord;
catch
  x_all_apical{:,i}=NaN;
y_all_apical{:,i}=NaN;  
end
end
%% Basal
for i=1:147
   
X_coord=[];
Y_coord=[];
try
[T X_coord Y_coord]=T_tree_sw (zz{:,i}{2}, '-s');
x_all_basal{:,i}=X_coord;
y_all_basal{:,i}=Y_coord;
catch
  x_all_basal{:,i}=NaN;
y_all_basal{:,i}=NaN;  
end
end
%% 
%% Basal
for i=1:147
   
X_coord=[];
Y_coord=[];
try
[T X_coord Y_coord]=T_tree_sw (zz{:,i}{1}, '-s');
x_all_apical{:,i}=X_coord;
y_all_apical{:,i}=Y_coord;
catch
  x_all_apical{:,i}=NaN;
y_all_apical{:,i}=NaN;  
end
end

%% Ocularity   
par=eye_sp(resp_id)
e1=find(par==0);
e2=find(par==1);
par2=(ang_inL4(resp_id,3)-ang_inL4(resp_id,1))*69
%par2=span(resp_id,3)-span(resp_id,6)

[statsout]=dual_barplot(par2,e1,e2,0);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})
%% 
  g1=[];
  g2=[];
   a=[];
 a=find(od_out_iviv(:,1)>0.0 & abs(od_out_iviv(:,3))>0 & r_sq>0.2); 
 par=od_out_iviv(a,2)
   figure;scatter(par,tot_inputL4(a,1),10,od_out_iviv(a,2),'filled')
   %% 
     g1=[];
  g2=[];
   a=[];
   cuty=0.5
g1=find(abs(od_out_iviv(:,3))<cuty)
g2=find(abs(od_out_iviv(:,3))>cuty)
 par=L5fr(:,1)

[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);hold on
%% 
a=find(r_sq>0.1 & abs(od_out_iviv(:,3))>0)
corr_plot(abs(tot_inputL4(a,1)-tot_inputL4(a,2)),abs(od_out_iviv(a,3)),od_out_iviv(a,3),{'EX','IN',''});ylabel('IN C_{x} (µm)','Color','b');xlabel('EX C_{x} (µm)','Color','r');
%% 

for i=1:147
%whole map
tmp1=ex_map_raw(3:end,:,i);
tmp2=in_map_raw(:,:,i);
% tmp1=ex_map(:,:,i);
% tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
%tot_input(i,:)=[sum(tmp1(:))/length(tmp1);sum(tmp2(:))/length(tmp2)];
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map_raw(3:5,:,i);
tmp2=in_map_raw(3:5,:,i);
%  tmp1=ex_map(3:5,:,i);
%  tmp2=in_map(3:5,:,i);
tot_inputL23(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map_raw(6:7,:,i);
tmp2=in_map_raw(6:7,:,i);
%  tmp1=ex_map(6:7,:,i);
%  tmp2=in_map(6:7,:,i);
tot_inputL4(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
end
%% 

%% Span barplot
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0); 
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
par=span(a,1)*69;
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})
%% 
par1=L5fr(a,2);
par2=od_out_iviv(a,1)
corr_plot(par1(s1),par2(s1),[],{'EX','IN',''})
