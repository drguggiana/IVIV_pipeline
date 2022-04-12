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
subplot(6,10,i);
hold on;plot(oris,str(resp_id(i)).TCpref,'--o');
y=str(resp_id(i)).TCpref;
ym=mean(y);
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
    if str(resp_id(i)).ODIpref>0;
    hold on;plot(str(resp_id(i)).fit_resp(1:360));
        yp=(str(resp_id(i)).fit_resp(1:360));
 yp_s=yp(oris+1);
    else str(resp_id(i)).ODIpref<0; 
        ;hold on;plot(str(resp_id(i)).fit_resp(361:end));
            yp=(str(resp_id(i)).fit_resp(361:end));
 yp_s=yp(oris+1);
    end
end
r2=1-(sum(sqrt((y-yp_s').^2))/sum(sqrt((y-ym).^2)));
r_square(i)=r2;
title([num2str(str(resp_id(i)).OSIpref) ' / ' num2str(r_square(i))]);
y=[];
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
%% Plotting example cells 143 and 122 in structure
%panel A
id_m=143;id_m2=122;
plot_morphologies(str,id_m,2,1,1,1);
plot_morphologies(str,id_m2,2,1,1,1);
%% %Plot the tuning curves and width from the folders
fit_ex1=[];fit_ex2=[];
%cell 143
load('D:\Postdoc_Margrie\Projects\L23\Paper1\Revision\example_peak_curves_cell143\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex1=Fit(1).ipsi.FittedDataOri
ind_tr1=peaks(1).mean_trial_ipsi
ind_tr1_sigma=peaks(1).mean_trial_ipsi;
%cell 122
load('D:\Postdoc_Margrie\Projects\L23\Paper1\Revision\example_peak_curves_cell121\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat')
fit_ex2=Fit(2).contra.FittedDataOri;
ind_tr2=peaks(2).mean_trial_contra;
ind_tr2_sigma=peaks(2).mean_trial_contra;
%get the baseline responsedor the fits 
[FittedData1, BaselineRsp1, PrefRsp1, PrefDir1, Sigma1, OppResp1, Error1, R21 ] = FIT_Carandini( str(id_m).TCpref );
[FittedData2, BaselineRsp2, PrefRsp2, PrefDir2, Sigma2, OppResp2, Error2, R22 ] = FIT_Carandini( str(id_m2).TCpref );
%% Plot panel B changes based on Reviewer 3
peak=max(fit_ex1);
find(fit_ex1==peak);
peak2=max(fit_ex2);
find(fit_ex2==peak2);
pori=[0:45:315];
ori=90-[1:1:180];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 200, 700, 250]);
subplot(1,3,1)
errorbar(pori,nanmean(ind_tr1,2),nanstd(ind_tr1,[],2)/sqrt(4),'-o','Color','k');hold on
hold on;box off
ylabel('dR/R0');xlabel('Direction (deg)');
xlim([0 315]);xticks([0:45:315]);xtickangle(45);
set(gca,'FontSize',10);
ylim([-5 60]);
text(10,60,['gOSI= ' num2str(round(od_out_iviv(id_m,1),2))]);

subplot(1,3,2)
hold on;box off
errorbar(pori,nanmean(ind_tr2,2),nanstd(ind_tr2,[],2)/sqrt(4),'-o','Color','b')
ylim([-5 60]);
text(10,60,['gOSI= ' num2str(round(od_out_iviv(id_m2,1),2))],'Color','b');ylabel('dR/R0')
xlim([0 315]);xticks([0:45:315]);xtickangle(45);xlabel('Direction (deg)');set(gca,'FontSize',10);

subplot(1,3,3)
plot(circshift(ori,find(fit_ex1==peak)-90),fit_ex1-BaselineRsp1,'k');
hold on;
plot(circshift(ori,find(fit_ex2==peak2)-90),fit_ex2-BaselineRsp2,'b');
hold on;box off;
xlabel('Distance from pref. orientation (deg)');
text(10,60,['TW= ' num2str(round(od_out_iviv(id_m,7),2))]);
hold on;box off
text(10,55,['TW= ' num2str(round(od_out_iviv(id_m2,7),2))],'Color','b');
ylim([-5 60]);
set(gca,'FontSize',10);
ylabel('dR/R0');
%% Read out morphology for all cells in structure and perform sholl analysis
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
m_res=[];
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1 & str(i).resp==1 & ~isnan(str(i).OSIpref)==1
        m_res(i)=1;     
    elseif ~isnan(str(i).morph(1))==1& str(i).iviv==1 & str(i).resp==0 
        m_res(i)=0;  
    else
       m_res(i)=NaN;
    end
end
morph_res_sub=[];
morph_res_sub=find(m_res==1);
%R2 criterion for good fitted cells
morph_res_sub2=morph_res_sub;
morph_res_sub2(find(r_sq(morph_res_sub2)<0.3))=[];
%% Plot the example sholl analyis for cell 143 and 122, PANEL F
temp=zz{143};
figure;
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(temp{1,1}, 20, '-s');
[sb, ddb, sdb, XP, YP, ZP, iD] = sholl_tree(temp{1,2}, 20, '-s');
temp=zz{122};
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
ylabel('Tuning Width','Color','k');set(gca,'FontSize',10)
%% TW branch points
corr_plot(morph_parameters(morph_res_sub2,4)',od_out_iviv(morph_res_sub2,7)',[],{'','',''});

%% gOSI
corr_plot(morph_parameters(morph_res_sub,4)',od_out_iviv(morph_res_sub,1)',[],{'','',''});
%% gOSI
corr_plot(max_s(morph_res_sub),od_out_iviv(morph_res_sub,1)',[],{'','',''});
%% TW
corr_plot(morph_parameters(morph_res_sub2,14)',od_out_iviv(morph_res_sub2,7)',[],{'','',''});
%% Check tuning width error based on Rev3

%DRAGO

% data=[];data=od_out_iviv(morph_res_sub2,7);
% pall=[];rall=[];
% psel=[];rsel=[];
% count=1
% for k=1:1000
% for i=1:length(data);
% data_m=[];tt1=[];tt2=[];
% data_m=data;
% ved=1:length(data);
% tt1=randperm(length(data),i);
% [tf, loc] = ismember(ved,tt1);
% tt2=find(tf==0);
% data_m(tt1)=data_m(tt1)+10;
% data_m(tt2)=data_m(tt2)-10;
% [r p]=corrcoef(max_s(morph_res_sub2)',data_m);
% psel(i)=p(2);
% rsel(i)=r(2);
% end
% pall(count,:)=psel;
% rall(count,:)=rsel;
% count=count+1;
% end
% figure;set(gcf,'color','w');
% subplot(2,1,1);histogram(rall(:));xlabel('Pearson Correlation R');
% subplot(2,1,2);histogram(pall(:));xlabel('p-Value');
%% %% Prepare corr matrix for S3
df=[morph_parameters(morph_res_sub,9)];
db=[morph_parameters(morph_res_sub,19)];
com=[];com=[morph_parameters(morph_res_sub,2) nanmax(df,[],2) dis_s(morph_res_sub)'  morph_parameters(morph_res_sub,4)  max_s(morph_res_sub)' ...
    morph_parameters(morph_res_sub,12) nanmax(db,[],2) dis_s_ba(morph_res_sub)'  morph_parameters(morph_res_sub,14) max_s_ba(morph_res_sub)'  pia_input(morph_res_sub') od_out_iviv(morph_res_sub,[1 2 3 6])]
G1=correlation_matrix(com,0);
%% subsample TW
df=[];db=[];
df=[morph_parameters(morph_res_sub2,9)];
db=[morph_parameters(morph_res_sub2,19)];
com=[];com=[morph_parameters(morph_res_sub2,2) nanmax(df,[],2) dis_s(morph_res_sub2)'  morph_parameters(morph_res_sub2,4)  max_s(morph_res_sub2)' ...
    morph_parameters(morph_res_sub2,12) nanmax(db,[],2) dis_s_ba(morph_res_sub2)'  morph_parameters(morph_res_sub2,14) max_s_ba(morph_res_sub2)'   pia_input(morph_res_sub2') od_out_iviv(morph_res_sub2,[7])]
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
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))' pia_input(morph_res_sub(a))]
for i=1:size(par_c,2)

    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),4)), par_c(:,i))
end
% Circular correlation for DRI
a=[]
a=find(od_out_iviv(morph_res_sub,2)>0.25 & r_sq(morph_res_sub)>0.3) ; 
par_c=[];
rho2=[];pval2=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))'  pia_input(morph_res_sub(a))]
for i=1:size(par_c,2)

    [rho2(i) pval2(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),5)), par_c(:,i))
end
% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% Plot correlation matrix for S3 panel C
%Gf=G1(11:14,1:10)
Gf=G1(12:15,1:11)
tG=[Gf;G2(12:end,1:11);c_cor_e;G1(15,1:11)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(tG(1:7,:));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]) 
xticklabels({'Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Horizontal extent','Dis peak branch','Nr. branch points','Peak number crossing', 'Pial depth'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'gOSI','gDSI','ODI','Ca_{peak}','TW','ORI','DIR','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);
%% /readout + plot circular correlation plots + sine fits
a=[];
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
%corr_plot(df(a),od_out_iviv(morph_res_sub(a),4),[],{'','',''});
% Extent APICAL vs ORI plus sine fit
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
% hold on;plot([106 106],[100 400],'--','Color',[0.7 0 0.4],'Linewidth',0.5);
% hold on;plot([53 53],[100 400],'--','Color',[0 0.5 0.5],'Linewidth',0.5);
hold on; p1 = [53 120]; p2 = [53 102];                         % Second Point
dp = p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',1,'Color',[0 0.5 0.5]);
hold on; p1 = [106 120]; p2 = [106 102];                         % Second Point
dp = p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',1,'Color',[0.7 0 0.4]);
scatter(od_out_iviv(morph_res_sub(a),4),df(a),7,'o','MarkerEdgeColor','k','MarkerFaceColor','k');xlabel('ORI (deg)');ylabel('Horizontal Extent (µm)');box off
xlim([0 180]);xticks(0:45:180);set(gca,'FontSize',10);text(60,400,['cc=' num2str(round(c_cor(1,2),3))]);
text(140,400,['p=' num2str(round(c_pva(1,2),3))]);set(gca,'FontSize',10);
% Extent BASAL vs ORI plus sine fit
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 250, 250]);
hold on; p1 = [53 120]; p2 = [53 102];                         % Second Point
dp = p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',1,'Color',[0 0.5 0.5]);
hold on; p1 = [106 120]; p2 = [106 102];                         % Second Point
dp = p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',1,'Color',[0.7 0 0.4]);
scatter(od_out_iviv(morph_res_sub(a),4),db(a),7,'o','MarkerEdgeColor','k','MarkerFaceColor','k');xlabel('ORI (deg)');ylabel('Horizontal Extent (µm)');box off
xlim([0 180]);xticks(0:45:180);set(gca,'FontSize',10);text(60,400,['cc=' num2str(round(c_cor(1,7),3))]);
text(140,400,['p=' num2str(round(c_pva(1,7),3))]);set(gca,'FontSize',10);
ylim([100 400]);
%% Compare using bins: REVIEWER
a=[];
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3 & ...
    abs(od_out_iviv(morph_res_sub,3))>0) ;
% binoori=reshape([str(morph_res_sub(a)).Ori],[2,length(morph_res_sub(a))])';
% binodir_delta=binoori(:,1)-binoori(:,2);
% binoOriA=binoori;
s1a=[];s2a=[];s1b=[];s2b=[]
%define sectors width
sector1=12; 
sector2=28;
%define midpoint of widows
midpoint1=106; 
midpoint2=53;
%define range
s1a=[midpoint1-sector1]
s1b=[midpoint1+sector2]
s2a=[midpoint2-sector2]
s2b=[midpoint2+sector1]

g1=[];g2=[];g3=[];g4=[];
%VS - find mono cells belonging to sectors 
% bi_limit=0.25;
% idx_mono=find(abs(od_out_iviv(morph_res_sub(a),3))>bi_limit)
% g1=find(od_out_iviv(a(idx_mono),5)>s1a & od_out_iviv(a(idx_mono),5)<s1b);
% g2=find(od_out_iviv(a(idx_mono),5)>s2a & od_out_iviv(a(idx_mono),5)<s2b);
% %g3=find(od_out_iviv(a(idx_mono),5)>s3a & od_out_iviv(a(idx_mono),5)<s3b);
% g3=find(od_out_iviv(a(idx_mono),5)>s3a | od_out_iviv(a(idx_mono),5)<s3b);
% %g4=find(od_out_iviv(a(idx_mono),5)>s4a & od_out_iviv(a(idx_mono),5)<s4b);
% g4=find(od_out_iviv(a(idx_mono),5)>s4a & od_out_iviv(a(idx_mono),5)<s4b);
% %VS - find bi cells belonging to sectors 
% idx_bi=find(abs(od_out_iviv(a,3))<bi_limit)
% tempIdx=find(binodir(idx_bi,1)>s1a & binodir(idx_bi,1)<s1b & binodir(idx_bi,2)>s1a & binodir(idx_bi,2)<s1b);
% g1=[idx_mono(g1); idx_bi(tempIdx)];
% tempIdx=find(binodir(idx_bi,1)>s2a & binodir(idx_bi,1)<s2b & binodir(idx_bi,2)>s2a & binodir(idx_bi,2)<s2b);
% g2=[idx_mono(g2); idx_bi(tempIdx)];
% tempIdx=find(binodir(idx_bi,1)>s3a & binodir(idx_bi,1)<s3b & binodir(idx_bi,2)>s3a & binodir(idx_bi,2)<s3b);
% g3=[idx_mono(g3); idx_bi(tempIdx)];
% tempIdx=find(binodir(idx_bi,1)>s4a & binodir(idx_bi,1)<s4b & binodir(idx_bi,2)>s4a & binodir(idx_bi,2)<s4b);
% g4=[idx_mono(g4); idx_bi(tempIdx)];
% par=[];
% %combine 315 and 135
% s1=[g1' g2'];
% %combine 45 and 225
% s2=[g3' g4'];
 g1=find(od_out_iviv(morph_res_sub(a),4)>s1a & od_out_iviv(morph_res_sub(a),4)<s1b);
 g2=find(od_out_iviv(morph_res_sub(a),4)>s2a & od_out_iviv(morph_res_sub(a),4)<s2b);

%UNCOMMENT FOR LAYER AND EX/IN, respectively 
par=df
%par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_exL4(a,3)-ang_exL4(a,1))*69;
%par=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%par=(ang_exL5(a,3)-ang_exL5(a,1))*69;
% par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
%SWITCH SIGN CAUSE OF FLIP

[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%ylim([-70 70]);yticks(-70:35:70);
%ylim([-150 200]);yticks(-150:100:200);

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
 %% Reviewer 1 overview
 a1=[];a2=[]; b=[];a=[];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3 & od_out_iviv(morph_res_sub,4)>18  & od_out_iviv(morph_res_sub,4)<68)
%a2=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3 & od_out_iviv(morph_res_sub,4)>150)
b=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3 & od_out_iviv(morph_res_sub,4)>71 & od_out_iviv(morph_res_sub,4)<121)
 %a=[a1;a2];
 color_id={[0.7 0 0.4],[0 0.5 0.5]};
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 500, 400]);
subplot(2,2,3)
for i=1:length(a)
  m=plot_tree(str(morph_res_sub(a(i))).morphtraces{1,1},[1 0 0],[0 str(morph_res_sub(a(i))).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{2};
         set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
end


subplot(2,2,1)
  for i=1:length(b)
  m=plot_tree(str(morph_res_sub(b(i))).morphtraces{1,1},[1 0 0],[0 str(morph_res_sub(b(i))).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{1};;
         set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
end
%  morph_parameters(morph_res_sub(a),9)
%  figure;
%   m=plot_tree(str(morph_res_sub(a(1))).morphtraces{1,1},[1 0 0],[0 str(morph_res_sub(a(1))).morphtraces{1,5} 0],[],1,'-b');hold on;
  
%Basael
subplot(2,2,4)
for i=1:length(a)
  m=plot_tree(str(morph_res_sub(a(i))).morphtraces{1,2},[1 0 0],[0 str(morph_res_sub(a(i))).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{2};
         set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 %axis off;
 grid on
 set(gcf,'color','w');
end
subplot(2,2,2)
  for i=1:length(b)
  m=plot_tree(str(morph_res_sub(b(i))).morphtraces{1,2},[1 0 0],[0 str(morph_res_sub(b(i))).morphtraces{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{1};
         set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 450]);
 axis off;
 set(gcf,'color','w');
  end

%% Sorted plot
 a=[];
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
exte=df(a);
 jj=[];
% [tt jj]=sort(df(a));
[tt jj]=sort(od_out_iviv(morph_res_sub(a),4),'ascend')
% oriu=od_out_iviv(morph_res_sub(a),4);
% oriu_s=oriu(jj);
%
hl=[];
exte_s=exte(jj)
hl=exte(jj)/max(exte(jj));
%hl=oriu(jj)/max(oriu(jj));
 
% cmap1='plasma';
% cmm=colormap(cmap1);
% cmm_don=downsample(cmm,8);
% cmap=cmm_don(jj,:)

for i=1:30
allc{i}=str(morph_res_sub(a(i))).morphtraces 
end
figure;
for i=1:30
subplot(1,30,i)
  m=plot_tree(allc{1,jj(i)}{1,1},[1 0 0],[0 allc{1,jj(i)}{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = [hl(i) 0 0];
         set(gca,'Ydir','reverse');
  xlim([-250 250]);
 ylim([0 300]);
 axis off;
 set(gcf,'color','w');
%title(num2str(oriu_s(i)))
%title(num2str(exte_s(i)))
end

% for i=1:30
% subplot(5,6,i)
%   m=plot_tree(str(morph_res_sub(a(i))).morphtraces{1,1},[1 0 0],[0 str(morph_res_sub(a(i))).morphtraces{1,5} 0],[],1,'-b');hold on;
%         m.EdgeColor = [1 0 0];
%          set(gca,'Ydir','reverse');
%   xlim([-350 350]);
%  ylim([0 450]);
%  %axis off;
%  set(gcf,'color','w');
% %title(num2str(df(i)))
% title(num2str(exte(i)))
% end
%% 
%df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=[];a=find(od_out_iviv(morph_res_sub,2)>0.1 & r_sq(morph_res_sub)>0.3) ; 
 jj=[];
% [tt jj]=sort(df(a));
[tt jj]=sort(od_out_iviv(morph_res_sub(a),2),'ascend')
exte=od_out_iviv(morph_res_sub(a),2);
hl=[];
hl=parula;
ak=[];ak=[1:16:256]
for i=1:length(jj)
allc{i}=str(morph_res_sub(a(i))).morphtraces 
end
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400,200]);
for i=13:16
    c=i+2;
subplot(1,4,i-12)
  m1=plot_tree(allc{1,jj(c)}{1,1},[1 0 0],[0 allc{1,jj(c)}{1,5} 0],[],1,'-b');hold on;
  m2=plot_tree(allc{1,jj(c)}{1,2},[1 0 0],[0 allc{1,jj(c)}{1,5} 0],[],1,'-b');hold on;
  m3=plot_tree(allc{1,jj(c)}{1,3},[1 0 0],[0 allc{1,jj(c)}{1,5} 0],[],1,'-b');hold on;
        m1.EdgeColor = hl(i+ak(i)-1,:)
        m2.EdgeColor = hl(i+ak(i)-1,:)
        m3.EdgeColor = hl(i+ak(i)-1,:)
         set(gca,'Ydir','reverse');
   xlim([-400 400]);
   ylim([-800 800]);
 axis off;
 set(gcf,'color','w');
%title(num2str(oriu_s(i)))
%title(num2str(exte_s(i)))
end



%% Example of elomgated and not elongated 
figure;
color_id={[0.7 0 0.4],[0 0.5 0.5]};
  m=plot_tree(allc{1,jj(1)}{1,1},[1 0 0],[0 allc{1,jj(1)}{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{2};
         set(gca,'Ydir','reverse');
         hold on
    m=plot_tree(allc{1,jj(25)}{1,1},[1 0 0],[0 allc{1,jj(15)}{1,5} 0],[],1,'-b');hold on;
        m.EdgeColor = color_id{1};
         set(gca,'Ydir','reverse');
         xlim([-250 250]);
 ylim([0 300]);
 axis off;
 set(gcf,'color','w');
%title(num2str(oriu_s(i)))
%title(num2str(exte_s(i)))



  %% Reviewer 2
 %1) TW
  close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 1000, 200]);
subplot(1,5,1);histogram(pia_input(morph_res_sub2),'BinLimits',[min(pia_input(morph_res_sub2)),max(pia_input(morph_res_sub2))],'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','Orientation','horizontal');
  ylabel('Pial depth (µm)');xlabel('Cell count');box off;

  subplot(1,5,3);scatter(max_s(morph_res_sub2),pia_input(morph_res_sub2),'ko','filled');
xlabel('Peak Nr. crossing','Color','k');ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10)
set(gcf,'color','w');xlim([0 15]);  title('Apical');
[r p]=corrcoef(max_s(morph_res_sub2),pia_input(morph_res_sub2))
text(2,300,['r=' num2str(round(r(2),3))]);
text(8,300,['p=' num2str(round(p(2),3))]);set(gca,'FontSize',10);

  subplot(1,5,2);scatter(morph_parameters(morph_res_sub2,4) ,pia_input(morph_res_sub2),'ko','filled');
xlabel('Nr. branch points','Color','k');ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10)
set(gcf,'color','w');  
 morph_parameters(morph_res_sub(a),4) 
 [r p]=corrcoef(morph_parameters(morph_res_sub2,4),pia_input(morph_res_sub2))
  text(4,300,['r=' num2str(round(r(2),3))]);
text(11,300,['p=' num2str(round(p(2),3))]);set(gca,'FontSize',10);
  
 subplot(1,5,4);
scatter(od_out_iviv(morph_res_sub2,1),pia_input(morph_res_sub2),'ko','filled');
xlabel('gOSI','Color','k');ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10)
set(gcf,'color','w');xlim([0 1])
[r p]=corrcoef(od_out_iviv(morph_res_sub2,1),pia_input(morph_res_sub2))
 text(0.2,300,['r=' num2str(round(r(2),3))]);
text(0.6,300,['p=' num2str(round(p(2),3))]);set(gca,'FontSize',10);


 subplot(1,5,5);
scatter(od_out_iviv(morph_res_sub2,7),pia_input(morph_res_sub2),'ko','filled');
xlabel('Tuning Width (deg)','Color','k');ylabel('Pial depth (µm)','Color','k');set(gca,'FontSize',10)
set(gcf,'color','w')
[r p]=corrcoef(od_out_iviv(morph_res_sub2,7),pia_input(morph_res_sub2))
xlim([0 40]);
 text(6,300,['r=' num2str(round(r(2),3))]);
text(22,300,['p=' num2str(round(p(2),3))]);set(gca,'FontSize',10);
%% Morphologu of unresponsive cells 

mg1=[];mg2=[];mo_unres=[];mo_res=[]
mg1=find(od_out_iviv(:,1)<0.25 & m_res'==1); 
mg2=find(m_res==0)';
%mo_unres=[mg1; mg2];
mo_unres=[mg2];
mo_res=find(od_out_iviv(:,1)>0.25 & m_res'==1); 

par=[];
par=max_s;
%par=max_s_ba;
%par=morph_parameters(:,4)

[statsout]=dual_barplot(par,mo_unres,mo_res,1);xticks([1:1:2]);hold on;xticklabels({'Unresponsive' ,'Tuned'});xtickangle(45);
ylabel('Peak nr. Sholl crossing');set(gca,'FontSize',10);
%ylim([-70 70]);yticks(-70:35:70);
title('Apical')

par=[];
par=max_s_ba;
%par=morph_parameters(:,4)

[statsout]=dual_barplot(par,mo_unres,mo_res,1);xticks([1:1:2]);hold on;xticklabels({'Unresponsive' ,'Tuned'});xtickangle(45);
ylabel('Peak nr. Sholl crossing');set(gca,'FontSize',10);
%ylim([-70 70]);yticks(-70:35:70);
title('Basal')
%% Sholl for aligned vs non_aligned
df=[];db=[];df=[morph_parameters(morph_res_sub,9)];db=[morph_parameters(morph_res_sub,19)];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ;
ty=[];ty=morph_res_sub(a);
%par=abs(ang_inL4(a,3)*69-ang_inL4(a,1)*69)
% par=ang_inL23(a,8)*69
% g1_1=find(od_out_iviv(ty,4)>0 & od_out_iviv(ty,4)<60) ;
% g1_2=find(od_out_iviv(ty,4)>150);
% g2_1=find(od_out_iviv(ty,4)>60 & od_out_iviv(ty,4)<150);
% g2=find(od_out_iviv(ty,4)>100 & od_out_iviv(ty,4)<150);
% [statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
% xticklabels({'10-60°','100-150°'});xtickangle(45);
% ylabel('CVL L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);



sa={};
dda={};
for i=1:length(ty)
    temp={}
    temp=zz{ty(i)};
    
[sa{i}, dda{i}, sd, XP, YP, ZP, iD] = sholl_tree(temp{1,1}, 20, '-s');
norm_sa{i}=sa{i}/max(sa{i});
end
%% 
colorMap = [linspace(0,1,256)', zeros(256,2)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 700, 450]);
Z=od_out_iviv(ty,4)/max(od_out_iviv(ty,4));
colorMap = [linspace(0,1,256)', zeros(256,2)]
for i=1:length(ty)
 hold on;
 %scatter(dda{i},norm_sa{i},25,repmat(od_out_iviv(ty(i),4),length(dda{i}),1),'filled')
 plot(dda{i},norm_sa{i},'color',[Z(i) 0 0],'LineWidth',1.5)
    
end
colormap(colorMap);
colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'0','45','90','135','180'});
     xlabel('Distance from Soma (microm)');
     ylabel('Normalized Nr. dendritic crossing');


    
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
%% 

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
set(gca,'XTick',[1.2536 4.8768 8.5 12.1232 15.7464],'XTickLabel',{'-500' '-250' '0'  '250' '500'});
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
set(gca,'XTick',[1.2536 4.8768 8.5 12.1232 15.7464],'XTickLabel',{'-500' '-250' '0'  '250' '500'});
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
set(gca,'XTick',[1.2536 4.8768 8.5 12.1232 15.7464],'XTickLabel',{'-500' '-250' '0'  '250' '500'});
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
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm) L2/3','Color','k');hold on;xlabel('Direction (deg)')
hold on;set(gca,'FontSize',10);ylim([-70 70]);yticks(-70:35:70);title('');xlim([0 360]);xticks(0:90:360);
%add estimated directions
hold on;plot([106 106],[-70 -50],'color',[0.7 0 0.4]);
hold on;plot([286 286],[-70 -50],'color',[0.7 0 0.4]);
hold on;plot([53 53],[-70 -50],'color',[0 0.5 0.5]);
hold on;plot([233 233],[-70 -50],'color',[0 0.5 0.5]);
hold on;text(25,70,'EX','color','r');
hold on;text(25,50,'IN','color','b');
%% REVIEWER
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Dir],[2,length(a)])';
binodir_delta=binodir(:,1)-binodir(:,2);
binoDirA=binodir;
% %define sectors width
% sector1=30; 
% sector2=50;
% %define midpoint of widows
% midpoint1=120; 
% midpoint2=30;
%define sectors width
sector1=24; 
sector2=56;
%define midpoint of widows
midpoint1=106; 
midpoint2=53;
%define range
s1a=[midpoint1-sector1]
s1b=[midpoint1+sector2]
s2a=[(midpoint1+180)-sector1]
s2b=[(midpoint1+180)+sector2]
s3a=357;
s3b=77;
s4a=[(midpoint2+180)-sector2]
s4b=[(midpoint2+180)+sector1]
g1=[];g2=[];g3=[];g4=[];
%VS - find mono cells belonging to sectors 
bi_limit=0.25;
idx_mono=find(abs(od_out_iviv(a,3))>bi_limit)
g1=find(od_out_iviv(a(idx_mono),5)>s1a & od_out_iviv(a(idx_mono),5)<s1b);
g2=find(od_out_iviv(a(idx_mono),5)>s2a & od_out_iviv(a(idx_mono),5)<s2b);
%g3=find(od_out_iviv(a(idx_mono),5)>s3a & od_out_iviv(a(idx_mono),5)<s3b);
g3=find(od_out_iviv(a(idx_mono),5)>s3a | od_out_iviv(a(idx_mono),5)<s3b);
%g4=find(od_out_iviv(a(idx_mono),5)>s4a & od_out_iviv(a(idx_mono),5)<s4b);
g4=find(od_out_iviv(a(idx_mono),5)>s4a & od_out_iviv(a(idx_mono),5)<s4b);
%VS - find bi cells belonging to sectors 
idx_bi=find(abs(od_out_iviv(a,3))<bi_limit)
tempIdx=find(binodir(idx_bi,1)>s1a & binodir(idx_bi,1)<s1b & binodir(idx_bi,2)>s1a & binodir(idx_bi,2)<s1b);
g1=[idx_mono(g1); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s2a & binodir(idx_bi,1)<s2b & binodir(idx_bi,2)>s2a & binodir(idx_bi,2)<s2b);
g2=[idx_mono(g2); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s3a | binodir(idx_bi,1)<s3b & binodir(idx_bi,2)>s3a | binodir(idx_bi,2)<s3b);
g3=[idx_mono(g3); idx_bi(tempIdx)];
tempIdx=find(binodir(idx_bi,1)>s4a & binodir(idx_bi,1)<s4b & binodir(idx_bi,2)>s4a & binodir(idx_bi,2)<s4b);
g4=[idx_mono(g4); idx_bi(tempIdx)];
par=[];
%combine 315 and 135
s1=[g1' g2'];
%combine 45 and 225
s2=[g3' g4'];
%UNCOMMENT FOR LAYER AND EX/IN, respectively 
%par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_exL4(a,3)-ang_exL4(a,1))*69;
par=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%par=(ang_exL5(a,3)-ang_exL5(a,1))*69;
 %par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
%SWITCH SIGN CAUSE OF FLIP
par(g2)=par(g2)*-1;
par(g4)=par(g4)*-1;
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%ylim([-70 70]);yticks(-70:35:70);
ylim([-150 200]);yticks(-150:100:200);

%% Absolute difference 
par=[];
%combine 315 and 135
s1=[g1' g2'];
%combine 45 and 225
s2=[g3' g4'];
%UNCOMMENT FOR LAYER AND EX/IN, respectively 
par=(ang_exL23(a,3)-ang_inL23(a,1))*69;
[statsout]=dual_barplot(abs(par),s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%% Centroid EX and IN 
 par1=[];par2=[];
   par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
   par2=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%     par1=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%   par2=(ang_exL4(a,3)-ang_exL4(a,1))*69;
%  par1=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%   par2=(ang_exL5(a,3)-ang_exL5(a,1))*69;
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
h=histogram(par2(s1)-par1(s1),3);box off;h.FaceColor=[0.5 0.5 0.5];
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
set(gca,'XTick',[1.2536 4.8768 8.5 12.1232 15.7464],'XTickLabel',{'-500' '-250' '0'  '250' '500'});



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
DSI_s2=od_out_iviv(a(s2),2);

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


[H,P]=ttest(rise,decay);


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


%% Alternative now test with gDSI as coulour after REVIEW
cmap1='plasma';
cm=colormap(cmap1);
cm_don=downsample(cm,round(length(cm)/length(DSI_s1)));
%colorMap = [linspace(0,1,length(DSI_s1))', zeros(length(DSI_s1),1) linspace(0,1,length(DSI_s1))'];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 200]);
data=[];data=riseDecayAligned_ex';
%Z=sort(DSI_s1/max(DSI_s1));
%Z=DSI_s1/max(DSI_s1)-0.1;
Z=[];
[ff, gg]=sort(DSI_s1)
Z=cm_don;
for i=1:length(Z)
%  hold on;
%  plot([1,2],[data(i,1),data(i,2)],'color',[Z(i) Z(i) Z(i)],'LineWidth',1) 
%  hold on; scatter([1,2],[data(i,1),data(i,2)],'MarkerFaceColor',[Z(i) Z(i) Z(i)],...
%      'MarkerEdgeColor',[Z(i) Z(i) Z(i)]) 
% hold on;pS=plotSpread([data(i,1),data(i,2)],'categoryIdx',[ones(1,length(data(i,1)))' ones(1,length(data(i,2)))'*2],...
%     'categoryMarkers',{'o','o'},'categoryColors',{colorMap(i,:),colorMap(i,:)});hold on;
 hold on;
 plot([1,2],[data(i,1),data(i,2)],'color',Z(i,:),'LineWidth',1) 
 hold on; scatter([1,2],[data(i,1),data(i,2)],'MarkerFaceColor',Z(i,:),...
     'MarkerEdgeColor',Z(i,:)) 
%  plot([1,2],[data(i,1),data(i,2)]) 
%   hold on; scatter([1,2],[data(i,1),data(i,2)]) 
end
xlim([0 3])
colormap(cm_don);
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);
xticks([1,2])
mm=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'0','0.25','0.5','0.75','1'});  set(get(mm,'title'),'string','gDSI');ylim([1 6.5])
     [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
%% Alternative now test with gDSI as coulour after REVIEW
%colorMap = [linspace(0,1,length(DSI_s2))', zeros(length(DSI_s2),1) linspace(0,1,length(DSI_s2))'];
cm_don=[];cmap1='plasma';
cm=colormap(cmap1);
cm_don=downsample(cm,round(length(cm)/length(DSI_s2)));
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 200]);
data=[];data=riseDecayNL';
Z=[];Z=cm_don;
for i=1:size(data,1)
 hold on;
 plot([1,2],[data(i,1),data(i,2)],'color',Z(i,:),'LineWidth',1) 
 hold on; scatter([1,2],[data(i,1),data(i,2)],'MarkerFaceColor',Z(i,:),...
     'MarkerEdgeColor',Z(i,:)) 
end
xlim([0 3])
colormap(cm_don);
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);
xticks([1,2])
mm=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'0','0.25','0.5','0.75','1'});set(get(mm,'title'),'string','gDSI');ylim([1 6.5])

%% alinged cells rise pref and null INHIBITION
cm_don=[];cmap1='plasma';
cm=colormap(cmap1);
cm_don=downsample(cm,round(length(cm)/length(DSI_s1)));
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 200]);
data=[];data=riseDecayAligned_in';
Z=[];Z=cm_don;
for i=1:size(data,1)
 hold on;
 plot([1,2],[data(i,1),data(i,2)],'color',Z(i,:),'LineWidth',1) 
 hold on; scatter([1,2],[data(i,1),data(i,2)],'MarkerFaceColor',Z(i,:),...
     'MarkerEdgeColor',Z(i,:)) 
end
xlim([0 3])
colormap(cm_don);
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);
xticks([1,2])
mm=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'0','0.25','0.5','0.75','1'});set(get(mm,'title'),'string','gDSI');ylim([1 6.5]);
     [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3)
%% Alternative now test with gDSI as coulour after REVIEW
%colorMap = [linspace(0,1,length(DSI_s2))', zeros(length(DSI_s2),1) linspace(0,1,length(DSI_s2))'];
cm_don=[];cmap1='plasma';
cm=colormap(cmap1);
cm_don=downsample(cm,round(length(cm)/length(DSI_s2)));
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 250, 200]);
data=[];data=riseDecayNL_in';
Z=[];Z=cm_don;
for i=1:size(data,1)
 hold on;
 plot([1,2],[data(i,1),data(i,2)],'color',Z(i,:),'LineWidth',1) 
 hold on; scatter([1,2],[data(i,1),data(i,2)],'MarkerFaceColor',Z(i,:),...
     'MarkerEdgeColor',Z(i,:)) 
end
xlim([0 3])
colormap(cm_don);
xticklabels({'Pref','Null'});ylabel('Length of rise');set(gca,'FontSize',10);
xticks([1,2])
mm=colorbar('Ticks',[0,0.25,0.5,0.75,1],...
         'TickLabels',{'0','0.25','0.5','0.75','1'});set(get(mm,'title'),'string','gDSI');ylim([1 6.5])

%% Correlation ratio of rise and DSI EX aligned
corr_plot([decayAligned_ex./riseAligned_ex]',[DSI_s1],[],{'','',''});ylim([0.25 1]);yticks([0.25:0.25:1]);
ylabel('gDSI');xlabel('Ratio rate of rise');set(gca,'Fontsize',10);xlim([0.4 1.8]);xticks([0.6:0.4:1.8]);
%% Correlation ratio of rise and DSI IN aligned
corr_plot([decayAligned_in./riseAligned_in]',[DSI_s1],[],{'','',''});ylim([0.25 1]);yticks([0.25:0.25:1]);
ylabel('gDSI');xlabel('Ratio rate of rise IN');set(gca,'Fontsize',10);xlim([0.4 1.8]);xticks([0.6:0.4:1.8]);
%% 

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
xlabel('EX-IN horizontal extent (µm)');
