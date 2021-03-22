%% START
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


%% FIGURE CREATOR



%% Figure 2: apical and basal complexity
%read ou morphology for all cells in structure and perform sholl analysis
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
%% Plot the example sholl analyis for cell 143 and 122, PANEL C
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
%% %Plot correlation tuning width vs apical and basal
corr_plot(max_s(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 15]);
ylabel('Tuning Width','Color','k');xlabel('Peak number of sholl crossings','Color','k');set(gca,'FontSize',10)
ylim([5 40]);yticks([0:10:40]);
%plot peak nr of sholl crossings BASAL vs TW, PANEL E
corr_plot(max_s_ba(morph_res_sub2)',od_out_iviv(morph_res_sub2,7),[],{'','',''});xlim([0 30]);
ylabel('Tuning Width','Color','k')
%% Supplement (S4)
%% Sholl analyos for 36 cells
close all;[max_s_sub dis_s_sub max_s_ba_sub dis_s_ba_sub]=sholl_analysis(zz,mcell,1);set(gcf,'color','w');box off;
legend('Apical','Basal');legend boxoff;ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%% Histogram morpholgy parameters iviv
stri={'Total Length (µm)','Horizontal extent', 'Dis. peak branch (µm)','Nr. branch points','Peak Nr. crossing',...
     'Total Length (µm)','Horizontal extent', 'Dis. peak branch (µm)','Nr. branch points','Peak Nr. crossing'}
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 200])
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
%ylabel('Counts');
 end
%% Prepare corr matrix
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
%% Plot correlation matrix
Gf=G1(11:14,1:10)
tG=[Gf;G2(11:end,1:10);c_cor_e;G1(14,1:10)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(tG(1:7,:));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]) 
xticklabels({'Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'gOSI','gDSI','ODI','Ca_{peak}','TW','ORI','DIR','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%% Rolling average setup/readout
a=[];
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
corr_plot(df(a),od_out_iviv(morph_res_sub(a),4),[],{'','',''});
%% extent vs ORI
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
scatter(df(a),od_out_iviv(morph_res_sub(a),4));ylabel('ORI (deg)');xlabel('Horizontal Extent (microm)');box off
ylim([0 180]);yticks(0:45:180)

%% Rolling average
m_res=[];
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1 
        m_res(i)=1;  
        hor_extent(i)=morph_parameters(i,9);
        hor_extent_ba(i)=morph_parameters(i,19);
    else
       m_res(i)=NaN;
       hor_extent(i)=NaN;
       hor_extent_ba(i)=NaN;
    end
end
%% Rolling average of horizontal elongation of apical vs basal
parameter_vector =  hor_extent;
parameter_vector2 = hor_extent_ba;
rolling_avg_display(str,parameter_vector2,parameter_vector2,45,1,0);hold on;
ylabel('Horizontal extent (µm)','Color','k');hold on;xlabel('Orientation preference (deg)')
hold on;set(gca,'FontSize',10);
%% 




%% Figure 3
%% Show vertical per layer for iviv cells
 [stats_g] = display_inputs_part2([frv(iviv_cells,:)],[frh(iviv_cells,:)],frv(iviv_cells,1:16)-frv(iviv_cells,17:end),frh(iviv_cells,1:16)-frh(iviv_cells,17:end),[]);
%%  Plot average iviv cells maps
plot_avg_maps(str,iviv_cells,ex_map,in_map,pia_input,10,0,[]);
%% Plot example maps
cnr=1;
plot_avg_maps(str,iviv_cells(cnr),ex_map,in_map,pia_input,1,0,[]);
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
 %% Group selective and unselective cells
g1=[];g2=[];
t1=find(od_out_iviv(:,1)<=0.25)
g2=find(od_out_iviv(:,1)>0.25)
g1=[unres_id t1'];
% par=L5fr(:,2);
% [statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);hold on;xticklabels({'UT' ,'TU'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% Plot maps from these groups
%%  Plot average iviv cells maps
plot_avg_maps(str,g1,ex_map_raw,in_map_raw,pia_input,10,0,[]);
plot_avg_maps(str,g2,ex_map_raw,in_map_raw,pia_input,10,0,[]);
%% Plot vertical profile selective vs unselective profiles
alpha1 = 0.3; 
alpha2 = 1; 
frac_v=frv;
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [300,100, 250, 300]);set(gcf,'color','w');
mexp=errorbar(nanmean(frac_v(g1,1:16)),1:16,nanstd(frac_v(g1,1:16))/sqrt(size(frac_v(g1,1:16),1)),'horizontal','Color',[1 0.3 0],'LineWidth',1.5);set(gca,'Ydir','reverse');
mexp.CapSize=3;

% Set transparency (undocumented)
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
 hold on;line([-1 1], [2.5 2.5],'Color','k','LineStyle','--');hold on;line([-1 1], [5.5 5.5],'Color','k','LineStyle','--');
 hold on;line([-1 1], [11 11],'Color','k','LineStyle','--'); hold on;line([-1 1], [7.5 7.5],'Color','k','LineStyle','--');
mexp2=errorbar(nanmean(frac_v(g1,17:end))*-1,1:16,nanstd(frac_v(g1,17:end))/sqrt(size(frac_v(g1,17:end),1)),'horizontal','Color',[0 0.9 1],'LineWidth',1.5);set(gca,'Ydir','reverse');
  
% Set transparency (undocumented)
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha1])
mexp2.CapSize=3;xlim([-0.4 0.4]);xticks([-0.4:0.2:0.4]);hold on;ylim([1 16]);yticks([1:5:16]);yticklabels({'1','6','11','16'});xlabel('Fraction of total input');set(gca,'FontSize',10);ylabel('Depth (µm)');set(gca,'FontSize',10)
hold on;

mexp=errorbar(nanmean(frac_v(g2,1:16)),1:16,nanstd(frac_v(g2,1:16))/sqrt(size(frac_v(g2,1:16),1)),'horizontal','Color',[1 0 0.2],'LineStyle',':','LineWidth',1.5);set(gca,'Ydir','reverse');
 alpha = 0.8;   
% Set transparency (undocumented)
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
 mexp.CapSize=3;
mexp2=errorbar(nanmean(frac_v(g2,17:end))*-1,1:16,nanstd(frac_v(g2,17:end))/sqrt(size(frac_v(g2,17:end),1)),'horizontal','Color',[0.3 0 1],'LineStyle',':','LineWidth',1.5);set(gca,'Ydir','reverse');
mexp2.CapSize=3;
alpha = 0.8;   
% Set transparency (undocumented)
set([mexp2.Bar, mexp2.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp2.Line.ColorData(1:3); 255*alpha2]);
box off;
yticklabels({'69','414','759','1104'});
%% Plot horizontal profile selective vs unselective profiles
alpha1 = 0.3; 
alpha2 = 1; 
frac_h=L23h;
frac_diffh=frac_h(:,1:16)-frac_h(:,17:end);
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 700, 350]);set(gcf,'color','w');
%EX and IN Horizontal
subplot(1,3,1);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color',[1 0.3 0],'LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color',[0 0.9 1],'LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color',[1 0 0.2],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color',[0.3 0 1],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
ylabel('Fraction of total input');set(gca,'FontSize',10);

frac_h=[];frac_h=L4h;
subplot(1,3,2);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color',[1 0.3 0],'LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color',[0 0.9 1],'LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color',[1 0 0.2],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color',[0.3 0 1],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);

frac_h=[];frac_h=L5h;
subplot(1,3,3);hold on;
mexp=errorbar(nanmean(frac_h(g1,1:16)),nanstd(frac_h(g1,1:16))/sqrt(size(g1,1)),'Color',[1 0.3 0],'LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])
hold on;mexp.CapSize=3;
hold on;mexp=errorbar(nanmean(frac_h(g1,17:end))*-1,nanstd(frac_h(g1,17:end))/sqrt(size(g1,1)),'Color',[0 0.9 1],'LineWidth',1.5); 
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha1])

hold on;mexp=errorbar(nanmean(frac_h(g2,1:16)),nanstd(frac_h(g2,1:16))/sqrt(size(g2,1)),'Color',[1 0 0.2],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(g2,17:end))*-1,nanstd(frac_h(g2,17:end))/sqrt(size(g2,1)),'Color',[0.3 0 1],'LineStyle',':','LineWidth',1.5);
set([mexp.Bar, mexp.Line], 'ColorType', 'truecoloralpha', 'ColorData', [mexp.Line.ColorData(1:3); 255*alpha2])
mexp.CapSize=3;ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');
xticklabels({'-552','-138','138','552'});
%% Difference horizontal ex and inh vs gOSI gDSI MEXICAN HAT THEORY
a=[];a=find(~isnan(od_out_iviv(:,1)));
corr_plot(span(a,1)-span(a,2),od_out_iviv(a,2),[],{'','',''});
%% Difference in fraction for ex and in gOSI, gDSI
a=[];a=find(~isnan(od_out_iviv(:,1)));
corr_plot(L23fr(a,1)-L23fr(a,2),od_out_iviv(a,1),[],{'','',''});


%% Figure 4
%% Rolling average
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm)','Color','k');hold on;xlabel('Direction preference (deg)')
hold on;set(gca,'FontSize',10)
%% Bar plot Cx for IN and EX in L23, L4, L5
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
%par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_exL4(a,3)-ang_exL4(a,1))*69;
%par=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%par=(ang_exL5(a,3)-ang_exL5(a,1))*69;
par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g2)=par(g2)*-1
par(s2)=par(s2)*-1
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
ylim([-70 70]);yticks(-70:35:70)
ylim([-150 200]);yticks(-150:50:200)
% ylim([-150 100]);yticks(-150:50:100)
%% Plot fraction for aligned and non aligned 
statsout=plot_horizontal_fraction_group(str,od_out_iviv,L23h,L4h,L5h)
%% Test ex and in differences
%% Centroid EX and IN 
% par1=[];par2=[];
%  par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%  par1(g2)=par1(g2)*-1;
%  par2=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%  par2(g2)=par2(g2)*-1;
par1=[];par2=[];
 par1=(ang_inL5(a,3)-ang_inL5(a,1))*69;
 par1(g2)=par1(g2)*-1;
 par2=(ang_exL5(a,3)-ang_exL5(a,1))*69;
 par2(g2)=par2(g2)*-1;
data=[par2(s1) par1(s1)];
p1=paired_plot([par2(s1) par1(s1)],0,{'r','b'});xticklabels({'EX','IN'});ylabel('C_{x} (µm)');set(gca,'FontSize',10)
[p1]=signrank(data(:,1) ,data(:,2),'tail','left')
%% histogram
fig7= figure;set(fig7, 'Name', 'Input distribution');set(fig7, 'Position', [200, 0, 200, 300]);set(gcf,'color','w');
h=histogram(par2(s1)-par1(s1),4);box off;h.FaceColor=[0.5 0.5 0.5];
xlabel('delta EX-IN');ylabel('Counts');set(gca,'FontSize',10);
hold on;line([0 0], [0 5],'Color','k','LineStyle','--');
hold on;scatter(nanmean(par2(s1)-par1(s1)),5,'kv','filled');
[h,p] = ttest(par2(s1)-par1(s1));
p = signrank(par2(s1)-par1(s1));
%% 
par3=od_out_iviv(a,2);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
scatter(par2(s1)-par1(s1),par3(s1));ylabel('DSI');xlabel('Cx EX - IN');box off
ylim([0 1]);yticks(0:0.25:1);set(gca,'FontSize',10);
%% 
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
%% 
par1=[];par2=[];
 par1=(ang_inL5(a,3)-ang_inL5(a,1))*69;
 par1(g2)=par1(g2)*-1;
 par2=(ang_exL5(a,3)-ang_exL5(a,1))*69;
 par2(g2)=par2(g2)*-1;
 %% 
 
par3=od_out_iviv(a,2);
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
scatter(par1(s1),par3(s1));ylabel('DSI');xlabel('Cx');box off
ylim([0 1]);yticks(0:0.25:1);set(gca,'FontSize',10);
