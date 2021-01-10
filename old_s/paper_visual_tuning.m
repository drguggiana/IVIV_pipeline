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
%% FIGURE 2
%% Show vertical per layer for iviv cells
 [stats_g] = display_inputs_part2([frv(iviv_cells,:)],[frh(iviv_cells,:)],frv(iviv_cells,1:16)-frv(iviv_cells,17:end),frh(iviv_cells,1:16)-frh(iviv_cells,17:end),[]);
%%  Plot average iviv cells maps
plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
%% Plot example maps
cnr=1;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);
%% Plot difference between vertical excitation and inhbition per cell
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
subplot(1,3,1);
histogram(diffL23fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL23fr),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L2/3');ylabel('Cell Counts');set(gca,'FontSize',10)
subplot(1,3,2);
histogram(diffL4fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL4fr),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L4');set(gca,'FontSize',10);xlabel('EX-IN vertical fraction');set(gca,'FontSize',10);
text(-0.5,28,'IN','Color','b');set(gca,'FontSize',10);text(0.5,28,'EX','Color','r');set(gca,'FontSize',10)
subplot(1,3,3);
histogram(diffL5fr(iviv_cells),'BinWidth',0.07,'FaceColor','w');box off;ylim([0 30]);xlim([-0.7 0.7]);xticks([-0.7:0.35:0.7]);yticks([0:15:30]);
hold on;scatter(nanmean(diffL5fr),30,'kv','filled')
hold on;line([0 0], [0 30],'Color','k','LineStyle','--');text(-0.5,24,'L5');set(gca,'FontSize',10);
%stats against 0
[p,h] = signrank(diffL5fr);
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
 %% Orientation tuning
 %% 
 
 
%% Direction tuning asymmetry
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
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);

%% Plot fraction for aligned and non aligned 
statsout=plot_horizontal_fraction_group(str,od_out_iviv,L23h,L4h,L5h)

