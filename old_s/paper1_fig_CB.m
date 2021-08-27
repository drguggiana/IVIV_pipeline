%Paper 1 figures for Current Biology 210812SW
%% 
%Load structure with 147 cells used for the paper
str      = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\old ivivstr';
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

% Raw input maps 
for i=1:length(str)
%whole map
tmp1=[];
tmp2=[];
tmp1=ex_map_raw(3:end,:,i);
tmp2=in_map_raw(:,:,i);
% tmp1=ex_map(:,:,i);
% tmp2=in_map(:,:,i);
tot_input(i,:)=[sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))];
%tot_input(i,:)=[sum(tmp1(:))/length(tmp1);sum(tmp2(:))/length(tmp2)];
%L1
tmp1=[];
tmp2=[];
tmp1=ex_map_raw(1:2,:,i);
tmp2=in_map_raw(1:2,:,i);
tot_inputL1raw(i,:)=abs([sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))]);
tmp1=[];
tmp2=[];
%L23
tmp1=ex_map_raw(3:5,:,i);
tmp2=in_map_raw(3:5,:,i);
tot_inputL23raw(i,:)=abs([sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))]);
tmp1=[];
tmp2=[];
%L4
tmp1=ex_map_raw(6:7,:,i);
tmp2=in_map_raw(6:7,:,i);
tot_inputL4raw(i,:)=abs([sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))]);
tmp1=[];
tmp2=[];
%L5
tmp1=ex_map_raw(8:11,:,i);
tmp2=in_map_raw(8:11,:,i);
tot_inputL5raw(i,:)=abs([sum(tmp1(:))/length(nonzeros(tmp1));sum(tmp2(:))/length(nonzeros(tmp2))]);
end

% Plot data with fits and calculate R2
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
% Add r2 square to str
for i=1:length(r_sq)
    str(i).r_sq=r_sq(i);
end

%  Extract horizontal fraction per layer
[L1h L23h L4h L5h] = h_fraclayer(ex_map, in_map);


%% %%%%%%%%%%%%%FIGURE CREATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 2: Apical but not basal dendritic morphology is related to orientation preference and tuning sharpness
%Panel A and B are currently produced by a different script, need to adapt
%in future, panel C is an entire scheme produced in Adobe Illustrator
%% Read ou morphology for all cells in structure and perform sholl analysis
close all;
for i=1:length(str)
    if ~isempty(str(i).morph)==1 
    zz{:,i}=str(i).morphtraces;
    else
        zz{:,i}=NaN;
    end
end
[max_s dis_s max_s_ba dis_s_ba]=sholl_analysis(zz,1:147,1);
%% Select all in vivo in vitro morpho cells (36)
mcell=find([str(:).iviv]==1 & morph_cells==1);
morph_sel=[];
%select 5 morpho paraemters
morph_sel=[morph_parameters(mcell,[2]) morph_parameters(mcell,[9]) dis_s(mcell)' morph_parameters(mcell,[4]) max_s(mcell)' ...
    morph_parameters(mcell,[12]) morph_parameters(mcell,[19]) dis_s_ba(mcell)' morph_parameters(mcell,[14]) max_s_ba(mcell)'];
%% Select all visual responsive cells
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1 & str(i).resp==1 & ~isnan(str(i).OSIpref)==1
        m_res(i)=1;     
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=find(m_res==1);
%R2 criterion for good fitted cells
morph_res_sub2=morph_res_sub;
morph_res_sub2(find(r_sq(morph_res_sub2)<0.3))=[];
%% Plot the example sholl analyis for cell 143 and 122, PANEL F
temp=zz{143};
figure;
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(temp{1,1}, 20, '-s');
[sb, ddb, sdb, XP, YP, ZP, iD] = sholl_tree(temp{1,2}, 20, '-s');
temp=zz{122}
[s1, dd1, sd1, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,1}, 20, '-s');
[s1b, dd1b, sd1b, XP1, YP1, ZP1, iD1] = sholl_tree(temp{1,2}, 20, '-s');
close(gcf);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 500, 250, 225]);
;plot(dd,s,'-k');box off;
hold on;plot(ddb,sb,'--k')
hold on;plot(dd1,s1,'-b');
hold on;plot(dd1b,s1b,'--b');
legend('Apical','Basal');legend boxoff;
ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%% %Plot correlation tuning width vs apical and basal PANEL G
corr_plot(max_s(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 15]);
ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%plot peak nr of sholl crossings BASAL vs TW, PANEL H
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 30]);
ylabel('Tuning Width','Color','k');
%% %% Prepare corr matrix for S3
df=[morph_parameters(morph_res_sub,9)];
db=[morph_parameters(morph_res_sub,19)];
com=[];com=[morph_parameters(morph_res_sub,2) nanmax(df,[],2) dis_s(morph_res_sub)'  morph_parameters(morph_res_sub,4)  max_s(morph_res_sub)' ...
    morph_parameters(morph_res_sub,12) nanmax(db,[],2) dis_s_ba(morph_res_sub)'  morph_parameters(morph_res_sub,14) max_s_ba(morph_res_sub)'  od_out_iviv(morph_res_sub,[1 2 3 6])]
G1=correlation_matrix(com,0);
%% subsample TW
df=[];db=[];
df=[morph_parameters(morph_res_sub2,9)];
db=[morph_parameters(morph_res_sub2,19)];
com=[];com=[morph_parameters(morph_res_sub2,2) nanmax(df,[],2) dis_s(morph_res_sub2)'  morph_parameters(morph_res_sub2,4)  max_s(morph_res_sub2)' ...
    morph_parameters(morph_res_sub2,12) nanmax(db,[],2) dis_s_ba(morph_res_sub2)'  morph_parameters(morph_res_sub2,14) max_s_ba(morph_res_sub2)'  od_out_iviv(morph_res_sub2,[7])]
G2=[];
G2=correlation_matrix(com,0);
%% Circular correlation for ORI
a=[];
df=[];db=[];
df=[morph_parameters(morph_res_sub,9)];
db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
par_c=[];rho1=[];pval1=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),4)), par_c(:,i))
end
% Circular correlation for DRI
a=[]
a=find(od_out_iviv(morph_res_sub,2)>0.25 & r_sq(morph_res_sub)>0.3) ; 
par_c=[];
rho2=[];pval2=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho2(i) pval2(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),5)), par_c(:,i))
end
% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% Plot correlation matrix for S3 panel C
Gf=G1(11:14,1:10)
tG=[Gf;G2(11:end,1:10);c_cor_e;G1(14,1:10)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(tG(1:7,:));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]) 
xticklabels({'Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'gOSI','gDSI','ODI','Ca_{peak}','TW','ORI','DIR','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%% Rolling average setup/readout + plot circular correlation plots + sine fits
a=[];
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
%corr_plot(df(a),od_out_iviv(morph_res_sub(a),4),[],{'','',''});
% Extent APICAL vs ORI plus sine fit
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
scatter(od_out_iviv(morph_res_sub(a),4),df(a),7,'o','MarkerEdgeColor','k','MarkerFaceColor','k');xlabel('ORI (deg)');ylabel('Horizontal Extent (µm)');box off
xlim([0 180]);xticks(0:45:180);set(gca,'FontSize',10);text(100,400,'cc=0.5');
text(150,400,'p<0.05');set(gca,'FontSize',10);
[m k]=sort(od_out_iviv(morph_res_sub(a),4))
tempn=df(a)
hold on;
[r2_ap x yp] = fit_sine(m',tempn(k)',1);
[mean_shuffle_ap CI_shuffle_ap] = shuffle_plot(m,tempn(k));
hold on;
text(25,400,['R2=' mat2str(round(r2_ap,2))],'Color','m');
set(gca,'FontSize',10);
% Extent BASAL vs ORI plus sine fit
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
scatter(od_out_iviv(morph_res_sub(a),4),db(a),7,'o','MarkerEdgeColor','k','MarkerFaceColor','k');xlabel('ORI (deg)');ylabel('Horizontal Extent (µm)');box off
xlim([0 180]);xticks(0:45:180);set(gca,'FontSize',10);text(100,400,'cc=0.19');
text(160,400,'n.s.');set(gca,'FontSize',10);
ylim([100 400]);
[m k]=sort(od_out_iviv(morph_res_sub(a),4))
tempn=db(a)
hold on
[r2_ba x yp] = fit_sine(m',tempn(k)',1);
[mean_shuffle_ba CI_shuffle_ba] = shuffle_plot(m,tempn(k));
hold on;
text(25,400,['R2=' mat2str(round(r2_ba,2))],'Color','m');
set(gca,'FontSize',10);
%% Compare to shuffle shown in S4 panel D
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 230, 250]);
b1=bar(1:2,[r2_ap r2_ba]); b1.FaceColor=[0.5 0.5 0.5];hold on;
errorbar(1:2,vertcat(mean_shuffle_ap,mean_shuffle_ba),vertcat(CI_shuffle_ap(1),CI_shuffle_ba(1)),vertcat(CI_shuffle_ap(2),CI_shuffle_ba(2)),'ok');
ylabel('R squared');xticklabels({'Apical','Basal'});box off;
set(gca,'FontSize',10);
%% Sholl analysis + plotting for 36 cells for S3 panel A 
close all;[max_s_sub dis_s_sub max_s_ba_sub dis_s_ba_sub]=sholl_analysis(zz,mcell,1);set(gcf,'color','w');box off;
legend('Apical','Basal');legend boxoff;ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%% Histogram morpholgy parameters iviv, plotting for S3 panel B
stri={'Total Length (µm)','Horizontal extent', 'Dis. peak branch (µm)','Nr. branch points','Peak Nr. crossing',...
     'Total Length (µm)','Horizontal extent', 'Dis. peak branch (µm)','Nr. branch points','Peak Nr. crossing'}
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200]);
 for i=1:10
hold on;
subplot(2,5,i)
h=histogram(morph_sel(:,i),8,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5)
h.EdgeColor = 'k';
h.FaceColor = [0.5 0.5 0.5];
xlabel(stri(i));
ylim([0 20]);
xlim([0 max(morph_sel(:,i))+max(morph_sel(:,i))*0.25]);
hAxis = gca;
hAxis.YAxisLocation = 'left';    % 'left' (default) or 'right'
hAxis.XAxisLocation = 'bottom'
box off
 end

 
 
%% Figure 3: Laminar and horizontal distributions of local excitatory and inhibitory inputs to functionally characterized cells are diverse
% Panel C and part of panel A is an entire scheme produced in Adobe Illustrator
%% Plot example maps for panel A
%cell 1, 35 and 70, use subfunction plot_avg_maps
cnr=1;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);
cnr=35;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);
cnr=70;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);
%% Show vertical per layer for iviv cells, use only top left plot for paper, panel B
 [stats_g] = display_inputs_part2([frv(iviv_cells,:)],[frh(iviv_cells,:)],frv(iviv_cells,1:16)-frv(iviv_cells,17:end),frh(iviv_cells,1:16)-frh(iviv_cells,17:end),[]);
%%  Plot average iviv cells maps, panel B
plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
%% Plot horizontal fraction per layer panel C
%using subfunction plot_horizontal_fraction
plot_horizontal_fraction(iviv_cells,L23h,L4h,L5h);
%% Plot horizontal centroid per layer with ref line and exclude outliers in L5, panel E
g=[];g=find(ang_inL4(iviv_cells,3)*69-ang_inL4(iviv_cells,1)*69<135);
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
 %% Group selective and unselective cells and plot vertical profile, panel F
g1=[];g2=[];
t1=find(od_out_iviv(:,1)<=0.25);
g2=find(od_out_iviv(:,1)>0.25);
g1=[unres_id t1'];g1=g1';
% Plot vertical profile selective vs unselective profiles
alpha1 = 0.7; 
alpha2 = 0.7; 
frac_v=frv;
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [300,100, 250, 300]);set(gcf,'color','w');
mexp=errorbar(nanmean(frac_v(g1,1:16)),1:16,nanstd(frac_v(g1,1:16))/sqrt(size(frac_v(g1,1:16),1)),'horizontal','Color','r','LineWidth',1.5);set(gca,'Ydir','reverse');
mexp.CapSize=3;
% Set transparency 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
 hold on;line([-1 1], [2.5 2.5],'Color','k','LineStyle','--');hold on;line([-1 1], [5.5 5.5],'Color','k','LineStyle','--');
 hold on;line([-1 1], [11 11],'Color','k','LineStyle','--'); hold on;line([-1 1], [7.5 7.5],'Color','k','LineStyle','--');
mexp2=errorbar(nanmean(frac_v(g1,17:end))*-1,1:16,nanstd(frac_v(g1,17:end))/sqrt(size(frac_v(g1,17:end),1)),'horizontal','Color','b','LineWidth',1.5);set(gca,'Ydir','reverse'); 
% Set transparency 
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha1])
mexp2.CapSize=3;xlim([-0.4 0.4]);xticks([-0.4:0.2:0.4]);hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});xlabel('Fraction of total input');set(gca,'FontSize',10);ylabel('Depth (µm)');set(gca,'FontSize',10)
hold on;
mexp=errorbar(nanmean(frac_v(g2,1:16)),1:16,nanstd(frac_v(g2,1:16))/sqrt(size(frac_v(g2,1:16),1)),'horizontal','Color','m','LineStyle',':','LineWidth',1.5);set(gca,'Ydir','reverse');
 alpha = 0.8;   
% Set transparency
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
 mexp.CapSize=3;
mexp2=errorbar(nanmean(frac_v(g2,17:end))*-1,1:16,nanstd(frac_v(g2,17:end))/sqrt(size(frac_v(g2,17:end),1)),'horizontal','Color','g','LineStyle',':','LineWidth',1.5);set(gca,'Ydir','reverse');
mexp2.CapSize=3;alpha = 0.8;   
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha2]);
box off;
yticklabels({'69','414','759','1104'});
% Stats for panel F
p=[];
for i=1:size(frac_v,2);
[p(i) k(i)]=ranksum(frac_v(g1,i),frac_v(g2,i));
end
nanmin(p)
%% %% Plot horizontal profile selective vs unselective profiles per layer, panel G
%note: display needs better coding, but works for now
alpha1 = 0.7; 
alpha2 = 0.7; 
%L23
L23hm=[];
frac_h=[];
L23hm=L23h;
frac_h=[];frac_h=L23hm;
[row,col] = find(frac_h<0);
for i=1:length(row)
    frac_h(row(i),col(i))=NaN;
    L23hm(row(i),col(i))=NaN;
end
frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color','r','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color','b','LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color','m','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color','g','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
ylabel('Fraction of total input');set(gca,'FontSize',10);

%%%%L4
L4hm=[];
L4hm=L4h;
frac_h=[];frac_h=L4hm;
[row,col] = find(frac_h<0);
for i=1:length(row)
    frac_h(row(i),col(i))=NaN;
    L4hm(row(i),col(i))=NaN;
end
subplot(1,3,2);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color','r','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color','b','LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color','m','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color','g','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);

%%%%%L5
L5hm=[];
L5hm=L5h;
frac_h=[];frac_h=L5hm;
[row,col] = find(frac_h<0);
for i=1:length(row)
    frac_h(row(i),col(i))=NaN;
    L5hm(row(i),col(i))=NaN;
end
subplot(1,3,3);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color','r','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color','b','LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color','m','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color','g','LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
% Statistics
p1=[];
for i=1:size(L23hm,2);
[p1(i) k]=ranksum(L23hm(g1,i),L23hm(g2,i)); 
end
p2=[];
for i=1:size(L4hm,2);
[p2(i) k]=ranksum(L4hm(g1,i),L4hm(g2,i)); 
end
p3=[];
for i=1:size(L5hm,2);
[p3(i) k]=ranksum(L5hm(g1,i),L5hm(g2,i)); 
end
%statistacally only compare more than 3 data points 
[sum(L23hm(g1,:)>0);sum(L23hm(g2,:)>0)];
[sum(L4hm(g1,:)>0);sum(L4hm(g2,:)>0)];
[sum(L5hm(g1,:)>0);sum(L5hm(g2,:)>0)];
nanmin([p1 p2 p3]);
%% Figure 4: Asymmetry of input distributions correlates with direction tuning 
% Panel A is an entire scheme produced in Adobe Illustrator
%% Plot of rolling average of direction vs Centroid Cx L23 
%subfunction rolling_avg_display needed
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm)','Color','k');hold on;xlabel('Direction preference (deg)')
hold on;set(gca,'FontSize',10);
%%  Barplot of centroid offset for alignend and nonaligned, panel C
a=[];s1=[]; s2=[];
%readout direction selective cells and cells with ODI greater than 0
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
binodir_delta=binodir(:,1)-binodir(:,2);
binoDirA=binodir;
%define sectors width
sector=80; 
%define midpoint of widows
midpoint=130; 
%define range
s1a=[midpoint-sector/2];s1b=[midpoint+sector/2];
s2a=[(midpoint+180)-sector/2];s2b=[(midpoint+180)+sector/2];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint+90)-sector/2];s4b=[(midpoint+90)+sector/2];
g1=[];g2=[];g3=[];g4=[];
%VS - find mono cells belonging to sectors 
bi_limit=0.25;
idx_mono=find(abs(od_out_iviv(a,3))>bi_limit)
g1=find(od_out_iviv(a(idx_mono),5)>s1a & od_out_iviv(a(idx_mono),5)<s1b);
g2=find(od_out_iviv(a(idx_mono),5)>s2a & od_out_iviv(a(idx_mono),5)<s2b);
g3=find(od_out_iviv(a(idx_mono),5)>s3a & od_out_iviv(a(idx_mono),5)<s3b);
g4=find(od_out_iviv(a(idx_mono),5)>s4a & od_out_iviv(a(idx_mono),5)<s4b);
%VS - find bi cells belonging to sectors 
idx_bi=find(abs(od_out_iviv(a,3))<bi_limit)
tempIdx=find(binodir(idx_bi,1)>s1a & binodir(idx_bi,1)<s1b & binodir(idx_bi,2)>s1a & binodir(idx_bi,2)<s1b);
g1=[idx_mono(g1); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s2a & binodir(idx_bi,1)<s2b & binodir(idx_bi,2)>s2a & binodir(idx_bi,2)<s2b);
g2=[idx_mono(g2); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s3a & binodir(idx_bi,1)<s3b & binodir(idx_bi,2)>s3a & binodir(idx_bi,2)<s3b);
g3=[idx_mono(g3); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s4a & binodir(idx_bi,1)<s4b & binodir(idx_bi,2)>s4a & binodir(idx_bi,2)<s4b);
g4=[idx_mono(g4); idx_bi(tempIdx)];
par=[];
%combine 315 and 135
s1=[g1' g2'];
%combine 45 and 225
s2=[g3' g4'];
%UNCOMMENT FOR LAYER AND EX/IN, respectively 
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_exL4(a,3)-ang_exL4(a,1))*69;
%par=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%par=(ang_exL5(a,3)-ang_exL5(a,1))*69;
% par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
%SWITCH SIGN CAUSE OF FLIP
par(g2)=par(g2)*-1;
par(g4)=par(g4)*-1;
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
ylim([-70 70]);yticks(-70:35:70);
%ylim([-150 200]);yticks(-150:100:200);

%% Centroid EX and IN 
 par1=[];par2=[];
%    par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%   par2=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%     par1=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%   par2=(ang_exL4(a,3)-ang_exL4(a,1))*69;
  par1=(ang_inL5(a,3)-ang_inL5(a,1))*69;
  par2=(ang_exL5(a,3)-ang_exL5(a,1))*69;
%  par2(g2)=par2(g2)*-1;
%par1=[];par2=[];
 %par1=(ang_inL5(a,3)-ang_inL5(a,1))*69;
 par1(g2)=par1(g2)*-1;
 %par2=(ang_exL5(a,3)-ang_exL5(a,1))*69;
 par2(g2)=par2(g2)*-1;
% data=[par2(s1) par1(s1)];
% p1=paired_plot([par2(s1) par1(s1)],0,{'r','b'});xticklabels({'EX','IN'});ylabel('C_{x} (µm)');set(gca,'FontSize',10)
% [p1]=signrank(data(:,1) ,data(:,2),'tail','left')
%% histogram
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 200, 300]);set(gcf,'color','w');
h=histogram(par2(s1)-par1(s1),4);box off;h.FaceColor=[0.5 0.5 0.5];
xlabel('delta EX-IN');ylabel('Counts');set(gca,'FontSize',10);
hold on;line([0 0], [0 5],'Color','k','LineStyle','--');
hold on;scatter(nanmean(par2(s1)-par1(s1)),5,'kv','filled');
[h,p] = ttest(par2(s1)-par1(s1));
p = signrank(par2(s1)-par1(s1))
%% Read out horizontal fraction per layer for sectors and plot profiles for panel D
frac_h=[];flip_hl=[];flip_h2=[];
frac_h=L23h(a,:);
flip_hl=[frac_h(g1,1:16); flip(frac_h(g2,1:16),2)];
flip_hl_in=[frac_h(g1,17:end); flip(frac_h(g2,17:end),2)];
flip_h2=[frac_h(g3,1:16); flip(frac_h(g4,1:16),2)];
flip_h2_in=[frac_h(g3,17:end); flip(frac_h(g4,17:end),2)];

%plot
%Excitation
fig3= figure;set(fig3, 'Name', 'comp hori');set(fig3, 'Position', [200, 300, 350, 400]);set(gcf,'color','w');hold on;
for i=1:size(flip_h2,1);
hold on;
p1=plot(flip_h2(i,1:16),'Color','k');
p1.Color(4)=0.25;
end
hold on;
for i=1:size(flip_hl,1);
hold on;
p1=plot(flip_hl(i,1:16),'r');
p1.Color(4)=0.25;
end
hold on;plot(nanmean(flip_h2),'Color','k','LineWidth',1.5);
hold on;plot(nanmean(flip_hl),'Color','r','LineWidth',1.5);box off;


%Inhibition 
hold on;
for i=1:size(flip_h2_in,1);
hold on;
p1=plot(flip_h2_in(i,1:16)*-1,'Color','k');
p1.Color(4)=0.25;
end
hold on;
for i=1:size(flip_hl_in,1);
hold on;
p1=plot(flip_hl_in(i,1:16)*-1,'b');
p1.Color(4)=0.25;
end
hold on;
plot(nanmean(flip_h2_in)*-1,'Color','k','LineWidth',1.5);
hold on;
plot(nanmean(flip_hl_in)*-1,'Color','b','LineWidth',1.5);box off;
hold on;line([8.5 8.5], [-0.35 0.35],'Color','k','LineStyle','--');
xlim([1 16]);ylim([-0.35 0.35]);yticks(-0.35:0.35:0.35);xticks(1:5:16);xticklabels({'-552','-138','138','552'});
xlabel('Horizontal position (µm)');ylabel('Fraction of total input');set(gca,'Fontsize',10);




%% Read out rate of rise etc for panel E, F and G

% Read out horizontal fraction per layer for sectors
frac_h=[];
%L23
frac_h=L23h(a,:);
%L4
% frac_h=L4h(a,:);
%L5
% frac_h=L5h(a,:);
%flip g1 which is 135 EX
hor_s1_ex=[flip(frac_h(g1,1:16),2); frac_h(g2,1:16)];
%sector 45/225 EX
hor_s2_ex=[frac_h(g3,1:16); flip(frac_h(g4,1:16),2)]; %VS19/3/21
%flip g1 which is 135 IN
hor_s1_in=[flip(frac_h(g1,17:32),2); frac_h(g2,17:32)];
%sector 45/225 IN
hor_s2_in=[frac_h(g3,17:32); flip(frac_h(g4,17:32),2)]; %VS19/3/21
% figure;plot(hor_s1_ex','--r');title('135 & 315');hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');hold on;plot(hor_s1_in'*-1,'--b');
% figure;plot(hor_s2_ex','--r');title('45 & 225');hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');hold on;plot(hor_s2_in'*-1,'--b');
DSI_s1=od_out_iviv(a(s1),2);

% rate of rise
rise=[];
decay=[];
riseIN=[];
decayIN=[];

for i=1:size(hor_s1_ex,1)
output=rise_decay_quantification(hor_s1_ex(i,:));
%output=rise_decay_quantification(hor_s1_in(i,:));
rise(i)=output.RHS_rise;
decay(i)=output.LHS_rise;
output=rise_decay_quantification(hor_s1_in(i,:));
riseIN(i)=output.RHS_rise;
decayIN(i)=output.LHS_rise;
end;
%EX
rise_mean_ex=mean(rise);
rise_std_ex=std(rise)/sqrt(length(rise));
decay_mean_ex=mean(decay);
decay_std_ex=std(decay)/sqrt(length(decay));
riseAligned_ex=rise;
decayAligned_ex=decay;
riseDecayAligned_ex=[rise; decay];

% rate of rise orthogonal
rise=[];
decay=[];
for i=1:size(hor_s2_ex,1)
output=rise_decay_quantification(hor_s2_ex(i,:));
rise(i)=output.RHS_rise;
decay(i)=output.LHS_rise;
end;

riseDecayNL=[rise; decay];
riseNL_mean=mean(rise);
riseNL_std=std(rise)/sqrt(length(rise));
decayNL_mean=mean(decay);
decayNL_std=std(decay)/sqrt(length(decay));
'rise decay orthogonal'
[H,P]=ttest(rise,decay)';


%IN
rise_mean_in=mean(riseIN);
rise_std_in=std(riseIN)/sqrt(length(riseIN));
decay_mean_in=mean(decayIN);
decay_std_in=std(decayIN)/sqrt(length(decayIN));
riseAligned_in=riseIN;
decayAligned_in=decayIN;
riseDecayAligned_in=[riseIN; decayIN];

% rate of rise orthogonal
rise=[];
decay=[];
riseIN=[];
decayIN=[];
for i=1:size(hor_s2_ex,1)
output=rise_decay_quantification(hor_s2_in(i,:));
rise(i)=output.RHS_rise;
decay(i)=output.LHS_rise;
end;

riseDecayNL_in=[rise; decay];
riseNL_mean_in=mean(rise);
riseNL_std_in=std(rise)/sqrt(length(rise));
decayNL_mean_in=mean(decay);
decayNL_std_in=std(decay)/sqrt(length(decay));

%% alinged cells rise pref and null EXCITATION
cl={'r','r'};
data=[];data=riseDecayAligned_ex';
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','r');    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','r','Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);

%% non alinged cells rise pref and null EXCITATION
cl={'k','k'};
data=[];data=riseDecayNL';
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color',[0.5 0.5 0.5]);    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','k','Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'NAL (M)','NAL (L)'});ylabel('Length of rise');set(gca,'FontSize',10);

%% alinged cells rise pref and null INHIBITION
cl={'b','b'};
data=[];data=riseDecayAligned_in';
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','b');    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','b','Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);

%% non alinged cells rise pref and null INHIBITION
cl={'k','k'};
data=[];data=riseDecayNL_in';
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color',[0.5 0.5 0.5]);    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','k','Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'NAL (M)','NAL (L)'});ylabel('Length of rise');set(gca,'FontSize',10);

%% Correlation ratio of rise and DSI EX aligned
corr_plot([decayAligned_ex./riseAligned_ex]',[DSI_s1],[],{'','',''});ylim([0.25 1]);yticks([0.25:0.25:1]);
ylabel('gDSI');xlabel('Ratio rate of rise EX');set(gca,'Fontsize',10);xlim([0.6 1.8]);xticks([0.6:0.4:1.8]);
%% Correlation ratio of rise and DSI IN aligned
corr_plot([decayAligned_in./riseAligned_in]',[DSI_s1],[],{'','',''});ylim([0.25 1]);yticks([0.25:0.25:1]);
ylabel('gDSI');xlabel('Ratio rate of rise IN');set(gca,'Fontsize',10);xlim([0.6 1.8]);xticks([0.6:0.4:1.8]);
