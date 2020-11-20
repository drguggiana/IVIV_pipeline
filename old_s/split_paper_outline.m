
%% START
%Load structure with 147 cells used for the paper
str_invitro       = 'D:\Postdoc_Margrie\Projects\L23\structure\';
folder_list = uipickfiles('FilterSpec',str_invitro);
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
 %% 
 
 
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
r2=1-(sum(sqrt((y-yp_s').^2))/sum(sqrt((y-ym).^2)))
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
%% 

close all;
% Histogram of pia distribution with color coded in vivo morpho by itself
% for panel B
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 300, 300]);
hold on;histogram(pia_input(iviv_cells),'FaceColor','k','EdgeColor','k','Orientation','horizontal');
hold on;histogram(pia_input(find([str(:).iviv]==1 & morph_cells==1)),'FaceColor','w','EdgeColor','k','Orientation','horizontal');
legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
    [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
ylim([100 350])
legend boxoff ;set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);
;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
%% 
%Example input maps of cell 99 in Panel D
plot_maps(str,ex_map_raw,in_map_raw,[1:147],99,pia_input);
%Polar plot peak normalized in Panel C
u=99;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [500, 600, 200, 200]);
polarplot(deg2rad(oris([1:end 1])),str(u).TuningCurve([1:8 1])/max(str(u).TuningCurve));hold on;
polarplot(deg2rad(oris([1:end 1])),str(u).TuningCurve([9:end 1])/max(str(u).TuningCurve));
ax = gca;rticks([0:0.5:1]);ax.LineWidth = 2;ax.ThetaDir = 'clockwise';ax.ThetaZeroLocation = 'left';
 ax.ThetaTick = [0:45:360]; 
 %% Plot average iviv cells maps
 plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
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
%% 

corr_plot((ang_exL23(iviv_cells,3)-ang_exL23(iviv_cells,1))*69,(ang_inL23(iviv_cells,3)-ang_inL23(iviv_cells,1))*69,[],{'EX','IN',''});ylabel('IN C_{x} (µm)','Color','b');xlabel('EX C_{x} (µm)','Color','r');

temp1=ang_exL4(iviv_cells,3)-ang_exL4(iviv_cells,1);
temp2=ang_inL4(iviv_cells,3)-ang_inL4(iviv_cells,1);
temp1(find(isnan(temp2)))=[];
temp2(find(isnan(temp2)))=[];

corr_plot(temp1*69,temp2*69,[],{'EX','IN',''});ylabel('IN C_{x} (µm)','Color','b');xlabel('EX C_{x} (µm)','Color','r');


%% 

temp1=[];temp2=[];
temp1=ang_exL5(iviv_cells,3)-ang_exL5(iviv_cells,1);
temp2=ang_inL5(iviv_cells,3)-ang_inL5(iviv_cells,1);
aa=[find(isnan(temp2)); find(isnan(temp1))]
temp1(aa)=[];
temp2(aa)=[];
corr_plot(temp1*69,temp2*69,[],{'EX','IN',''});ylabel('IN','Color','b');xlabel('EX','Color','r');xlim([-200 200]);ylim([-200 200]);xlabel('EX C_{x} (µm)','Color','r');
%% 

g=find(ang_inL4(iviv_cells,3)*69-ang_inL4(iviv_cells,1)*69<135)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
plot(ang_exL23(iviv_cells,3)*69-ang_exL23(iviv_cells,1)*69,ang_inL23(iviv_cells,3)*69-ang_inL23(iviv_cells,1)*69,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-100 100];
ref.YData=[-100 100];
ref.Color='k';box off;xlim([-115 100]);ylim([-115 100]);hold on;title('L23','FontWeight','normal');xticks([-100:50:100]);yticks([-100:50:100]);
subplot(1,3,2);
plot(ang_exL4(iviv_cells,3)*69-ang_exL4(iviv_cells,1)*69,ang_inL4(iviv_cells,3)*69-ang_inL4(iviv_cells,1)*69,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'; ref.XData=[-150 150];
ref.YData=[-150 150];ylabel('IN C_{x} (µm)','Color','b')
ref.Color='k';box off;xlim([-175 150]);ylim([-175 150]);hold on;title('L4','FontWeight','normal');xticks([-150:75:200]);yticks([-150:75:150]);
xlabel('EX C_{x} (µm)','Color','r')
subplot(1,3,3);
plot(ang_exL5(iviv_cells,3)*69-ang_exL5(iviv_cells,1)*69,ang_inL5(iviv_cells,3)*69-ang_inL5(iviv_cells,1)*69,'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);
set(gcf,'color','w');ref= refline(1,0);set(gca,'FontSize',10);ref.LineStyle='--'
ref.Color='k';box off;xlim([-335 300]);ylim([-335 300]);hold on;title('L5','FontWeight','normal');xticks([-300:150:300]);yticks([-300:150:300]);

ref.XData=[-300 300];
ref.YData=[-300 300];
%% Rolling average
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm)','Color','k');hold on;xlabel('Direction preference (deg)')
hold on;set(gca,'FontSize',10)
%% 
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
figure;scatter(abs((ang_inL23(a,3)-ang_inL23(a,1))*69),od_out_iviv(a,5))

[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_inL23(a,3)-ang_inL23(a,1))*69))
%% 

corr_plot(od_out_iviv(a,4)-45,abs((ang_inL23(a,3)-ang_inL23(a,1))*69),abs(od_out_iviv(a,3)),{'Orientation along slice','IN',''});
%% 
close all
fig3= figure;set(fig3, 'Name', 'Input distribution');set(fig3, 'Position', [200, 600, 200, 200]);set(gcf,'color','w');
scatter(od_out_iviv(a,4)-45,abs((ang_inL23(a,3)-ang_inL23(a,1))*69),12,'bo','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');

[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_inL23(a,3)-ang_inL23(a,1))*69));
text(-20,90,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','b')

hold on
scatter(od_out_iviv(a,4)-45,abs((ang_exL23(a,3)-ang_exL23(a,1))*69),12,'ro','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');

[rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs((ang_exL23(a,3)-ang_exL23(a,1))*69));
text(-20,100,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','r')
set(gca,'FontSize',10);
ylim([-10 100])
%% 

a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
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
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% Correlation with DSI

par2=[];
par2=od_out_iviv(a,2);
corr_plot(abs(par(s1)),par2(s1),[],{'PC1ex','PC1com','Pial depth'});
%% Correlation comparison across layers
a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
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
%% 
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
par=[];
par=L23h(a,10);
dd=[par(s1);par(s2)]
gr=[ones(length(par(s1)),1);ones(length(par(s2)),1)*2]
figure;
comb=[dd gr];
for i=1:size(comb,2)-1   
[p,tbl,stats] = kruskalwallis(comb(:,i),comb(:,end),'off');
multicom(:,:,i) = multcompare(stats,'CType','bonferroni');
signi_ev(i)=multicom(:,6,i)<0.05;
end

%% Plot horizontal fraction per layer
frac_h=[];
frac_h=L23h(a,:);
frac_diffh=[];
frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 600, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
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
frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
%EX and IN Horizontal
subplot(1,3,2);hold on;
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
frac_diffh=frac_h(s1,1:16)-frac_h(s1,17:end);
%EX and IN Horizontal
subplot(1,3,3);hold on;
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
%% 
miviv=find([str(:).iviv]==1 & morph_cells==1);
smo=morph_cells(a)
par=[];
par=morph_parameters(a,21)
[statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% 
%% Morphology

a=[];
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
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

s1=[g1' g2'];
s2=[g3' g4'];

par1=ap_map(:,:,a);

s1t=[];
par1(:,:,g1)=flip(par1(:,:,g1));
s1t=nansum(nansum(par1(:,1:8,s1)))
s1M=s1t(:);
s1t=[];
s1t=nansum(nansum(par1(:,9:end,s1)))
s1L=s1t(:);
[s1M s1L]
p1=paired_plot([s1M s1L],1,{'k','r'});xticklabels({'M','L'});ylabel('Density apical');;set(gca,'FontSize',10)

%% Morphology


par2=ba_map(:,:,a);
s1t=[];
par1(:,:,g1)=flip(par2(:,:,g1));
s1t=nansum(nansum(par2(:,1:8,s1)))
s1M=s1t(:);
s1t=[];
s1t=nansum(nansum(par2(:,9:end,s1)))
s1L=s1t(:);
[s1M s1L]

p1=paired_plot([s1M s1L],1,{'k','r'});xticklabels({'M','L'});ylabel('Density basal');;set(gca,'FontSize',10)
%% 
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

%% 
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
%% 

a=[];
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0); 
sector=60;
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
par=dis_s_ba(a);
[statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'NAL','AL'})