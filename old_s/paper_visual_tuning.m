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
%% Raw input maps 
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
% close all;fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [100, 600, 300, 300]);
% hold on;histogram(pia_input(iviv_cells),'FaceColor','k','EdgeColor','k','Orientation','horizontal');
% hold on;histogram(pia_input(find([str(:).iviv]==1 & morph_cells==1)),'FaceColor','w','EdgeColor','k','Orientation','horizontal');
% legend([' in vivo + input  (n=' num2str(length(pia_input(iviv_cells))),')'],...
%     [' Morphology (n=' num2str(length(pia_input(find([str(:).iviv]==1 & morph_cells==1)))),')']);
% ylim([100 350])
% legend boxoff ;set(gca,'Ydir','reverse');yticks([100:100:400]);set(gca,'FontSize',10);
% ;ylabel('Pial depth (µm)');xlabel('Cell count');box off;
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
%% Raw input and selectivity 
close all;
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 900, 200]);
subplot(1,4,1);
scatter(tot_inputL1raw(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;text(0.4,1,'L1');set(gca,'FontSize',10);xlabel('PSC amplitude (pA)');xlim([-5 round(max(max(tot_inputL1raw)),-1)]);ylabel('gOSI');
subplot(1,4,2);
scatter(tot_inputL23raw(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL23raw(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);
hold on;text(0.4,1,'L2/3');set(gca,'FontSize',10);xlabel('PSC amplitude (pA)');xlim([-5 round(max(max(tot_inputL23raw)),-1)])
subplot(1,4,3);
scatter(tot_inputL4raw(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL4raw(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');ylim([0 1]);xlabel('PSC amplitude (pA)');
hold on;text(0.3,1,'L4');set(gca,'FontSize',10);%xlim([-0.02 0.6])
xlim([-5 round(max(max(tot_inputL4raw)),-1)])
subplot(1,4,4);
scatter(tot_inputL5raw(:,2),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','b');
hold on;scatter(tot_inputL5raw(:,1),od_out_iviv(:,1),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','r');xlabel('PSC amplitude (pA)');ylim([0 1]);%xlim([-0.02 0.6])
hold on;text(0.3,1,'L5');set(gca,'FontSize',10);legend('IN','EX');
legend boxoff;
xlim([-5 round(max(max(tot_inputL5raw)),-1)])
%% 


%% Group selective and unselective cells
g1=[];g2=[];
t1=find(od_out_iviv(:,1)<=0.25)
g2=find(od_out_iviv(:,1)>0.25)
g1=[unres_id t1'];
par=L5fr(:,2);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);hold on;xticklabels({'UT' ,'TU'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
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

%% 


hold on;mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
xlabel('Horizontal map position')
hold on;text(2.5,0.6,'L2/3','FontSize',12);
frac_h=L4h(iviv_cells,:);hold on;
%% 



%% Direction tuning asymmetry

%% Rolling average
parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);hold on;ylabel('C_{x} (µm)','Color','k');hold on;xlabel('Direction preference (deg)')
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
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
%par=(ang_inL23(a,3)-ang_inL23(a,1))*69;
%par=(ang_exL4(a,3)-ang_exL4(a,1))*69;
%par=(ang_inL4(a,3)-ang_inL4(a,1))*69;
%par=(ang_exL5(a,3)-ang_exL5(a,1))*69;
%par=(ang_inL5(a,3)-ang_inL5(a,1))*69;
%par=od_out_iviv(a,3)
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL' ,'NAL'});ylabel('C_{x} (µm)');set(gca,'FontSize',10);xtickangle(45);
%ylim([-70 70]);yticks(-70:35:70)
ylim([-150 200]);yticks(-150:50:200)
%% Getting the cell which are orientation selective but non direction slective
b=[];
%OSI>0.25 and DSI<0.25
b=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)<0.25 & abs(od_out_iviv(:,3))>0); 
%remove binocular cells that are oppsoite
binodir=[];binodir_delta=[];
binodir=reshape([str(b).Ori],[2,length(b)])';
binodir_delta=binodir(:,1)-binodir(:,2);
idx_bi=[]; idx_ov=[];
idx_bi=find(abs(od_out_iviv(b,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(b,3))<0.1)))>45);
b(idx_bi(id_ov))=[];
%chose sctors for orientation
sector=45;
midpoint=130;
s3a=[];s4a=[];s3b=[];s4b=[];
s3a=[(midpoint-90)-sector/2];s3b=[(midpoint-90)+sector/2];
s4a=[(midpoint)-sector/2];s4b=[(midpoint)+sector/2];
g5=[];
g6=[];
%g5 contain the desired cells 
g5=find(od_out_iviv(b,4)>s3a & od_out_iviv(b,4)<s3b);
g6=find(od_out_iviv(b,4)>s4a & od_out_iviv(b,4)<s4b);
%
%% Combine these cells eith the prvious DSI>0.25 ones
%INDICES ARE FO THE 147 now again
idx_dsi_align=[a(s1); b(g5)];

%% Centroid EX and IN 
 par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
 par1(g1)=par1(g1)*-1;
 par2=(ang_exL23(a,3)-ang_exL23(a,1))*69;
 par2(g1)=par2(g1)*-1;
% par1=(ang_inL4(a,3)-ang_inL4(a,1))*69;
% par1(g1)=par1(g1)*-1;
% par2=(ang_exL4(a,3)-ang_exL4(a,1))*69;
% par2(g1)=par2(g1)*-1;
data=[par2(s1) par1(s1)];
p1=paired_plot([par2(s1) par1(s1)],0,{'r','b'});xticklabels({'EX','IN'});ylabel('C_{x} (µm)');set(gca,'FontSize',10)
[p1]=signrank(data(:,1) ,data(:,2),'tail','left')
%% Plot fraction for aligned and non aligned 
statsout=plot_horizontal_fraction_group(str,od_out_iviv,L23h,L4h,L5h)




%% Morphology horizontal sholl analysis from VS
%align all morhologies correctly using m_flip again (morph_flip is just a
%147 matlab table with 1, 0 and NaN)
morph_flip     = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\old ivivstr';
folder_list = uipickfiles('FilterSpec',morph_flip);
load(char(folder_list));
%% load cyl_coord from Dropbox
cd('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Manuscript\split_paper\Paper 1\Volker_xy_coord_branchpoint');
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
%% Plot exmaple cell for horizontal sholl analysis
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
%% Read out iviv cells with morph
for i=1:length(str)
    if ~isnan(str(i).morph(1))==1  & str(i).iviv==1
        m_iviv(i)=1;     
    else
       m_iviv(i)=NaN;
    end
end
m_iviv_idx=find(m_iviv==1);
%% Plot all cells with horizontal sholl analysis

figure;hold on;
for cr=1:length(m_iviv_idx)
%     scatter(mx_basal{1, m_iviv_idx(cr)}(:,1),mx_basal{1, m_iviv_idx(cr)}(:,3));
% hold on;scatter(mx_apical{1, m_iviv_idx(cr)}(:,1),mx_apical{1, m_iviv_idx(cr)}(:,3));
% set(gcf,'color','w');box off;
% axis off;set(gca, 'YDir','reverse');
  data=mx_apical{1, m_iviv_idx(cr)}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
 hold on; plot(output.resampleX,output.resampleScholl,'Color',[0.5 0.5 0.5]);
data=mx_basal{1, m_iviv_idx(cr)}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
plot(output.resampleX,output.resampleScholl,'m');
  set(gcf,'color','w');box off;ylabel('Horizontal sholl: Branch points');

end
legend({'Apical','Basal'});legend boxoff
%% Same for tuned and untuned

g1=[];g2=[];
t1=find(od_out_iviv(:,1)<=0.25) 
g2=find(od_out_iviv(:,1)>0.25)
g1=[unres_id t1'];
mi_ut=intersect(g1,m_iviv_idx);
mi_t=intersect(g2,m_iviv_idx);

figure;hold on;
hold on;
for cr=1:length(mi_t)
   data=mx_apical{1, mi_t(cr)}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
 hold on; plot(output.resampleX,output.resampleScholl,'Color','m');
data=mx_basal{1, mi_t(cr)}(:,1:2);

end
hold on;
for cr=1:length(mi_ut)
%     scatter(mx_basal{1, m_iviv_idx(cr)}(:,1),mx_basal{1, m_iviv_idx(cr)}(:,3));
% hold on;scatter(mx_apical{1, m_iviv_idx(cr)}(:,1),mx_apical{1, m_iviv_idx(cr)}(:,3));
% set(gcf,'color','w');box off;
% axis off;set(gca, 'YDir','reverse');
  data=mx_apical{1, mi_ut(cr)}(:,1:2);
  output=scholl_analysis_horizontal_vs1(data,5,300);
 hold on; plot(output.resampleX,output.resampleScholl,'Color','g');
data=mx_basal{1, mi_ut(cr)}(:,1:2);
  set(gcf,'color','w');box off;ylabel('Horizontal sholl: Branch points');
end

legend({'Tuned','Untuned'});legend boxoff
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
%% Plot aligned and non aligned groups APICAL
s1_out_dira=shell_scholl_analysis_horizontal_vs1(s1_coord,5,300,1);
s2_out_dira=shell_scholl_analysis_horizontal_vs1(s2_coord,5,300,1);
close all;
%% %% Plot aligned and non aligned groups BASAL
s1_out_dirb=shell_scholl_analysis_horizontal_vs1(s1_coord,5,300,2);
s2_out_dirb=shell_scholl_analysis_horizontal_vs1(s2_coord,5,300,2);
close all;
%% %% Plot averages with shaded errorbar 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 250]);
subplot(1,2,2);
shadedErrorBar(s2_out_dirb.resampleX,s2_out_dirb.meanScholl,s2_out_dirb.stdScholl,'lineProps','g');
hold on;
shadedErrorBar(s1_out_dirb.resampleX,s1_out_dirb.meanScholl,s1_out_dirb.stdScholl,'lineProps','m');
legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Basal tree');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10)
legend('NonAligned','Aligned');xlabel('Distance from Soma (µm)');
subplot(1,2,1);
shadedErrorBar(s2_out_dira.resampleX,s2_out_dira.meanScholl,s2_out_dira.stdScholl,'lineProps','g');
hold on;
shadedErrorBar(s1_out_dira.resampleX,s1_out_dira.meanScholl,s1_out_dira.stdScholl,'lineProps','m');
%legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Apical tree');
ylabel('Horizontal sholl: Branch points');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10)


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
%% Plot aligned and non aligned groups basal
s1_out_orb=shell_scholl_analysis_horizontal_vs1(s1_coord_or,5,300,2)
s2_out_orb=shell_scholl_analysis_horizontal_vs1(s2_coord_or,5,300,2);
close all;
%% Plot aligned and non aligned groups apical
s1_out_ora=shell_scholl_analysis_horizontal_vs1(s1_coord_or,5,300,1)
s2_out_ora=shell_scholl_analysis_horizontal_vs1(s2_coord_or,5,300,1);
close all;
%% Plot averages with shaded errorbar 
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 700, 250]);
subplot(1,2,2);
shadedErrorBar(s2_out_orb.resampleX,s2_out_orb.meanScholl,s2_out_orb.stdScholl,'lineProps','m');
hold on;
shadedErrorBar(s1_out_orb.resampleX,s1_out_orb.meanScholl,s1_out_orb.stdScholl,'lineProps','g');
legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Basal tree');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10)
legend('Aligned (135 deg)','NonAligned (45 deg)');xlabel('Distance from Soma (µm)');
subplot(1,2,1);
shadedErrorBar(s2_out_ora.resampleX,s2_out_ora.meanScholl,s2_out_ora.stdScholl,'lineProps','m');
hold on;
shadedErrorBar(s1_out_ora.resampleX,s1_out_ora.meanScholl,s1_out_ora.stdScholl,'lineProps','g');
%legend('AL','NAL');legend boxoff;
set(gcf,'color','w');box off;
title('Apical tree');
ylabel('Horizontal sholl: Branch points');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10)


%% Statistics using all complexity points from -80 to +80 basal
% sect=find(s2_out_or.resampleX>20 & s2_out_or.resampleX<80);  
% alls1b=s1_out_orb.tempScholl(:,sect);
% alls2b=s2_out_orb.tempScholl(:,sect);
% % [p r]=ranksum(alls1b(:),alls2b(:));
% %% Statistics using all complexity points from -80 to +80 apical
% sect=find(s2_out_ora.resampleX>20 & s2_out_ora.resampleX<80);  
% alls1a=s1_out_ora.tempScholl(:,sect);
% alls2a=s2_out_ora.tempScholl(:,sect);
% % [p2 r2]=ranksum(alls1a(:),alls2a(:))

%% Check whether different 
p = ranksum(s1_out_ora.tempScholl(:),s2_out_ora.tempScholl(:))
p = ranksum(s1_out_orb.tempScholl(:),s2_out_orb.tempScholl(:))

%% 

dg=[max(alls2a,[],2); max(alls1a,[],2)]
gg=[ones(1,length(max(alls2b,[],2)))'; ones(1,length(max(alls1b,[],2)))'*2] 
s1=find(gg==1) 
s2=find(gg==2) 
par=dg;
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL','NAL'});ylabel('Max horizontal extent apical')
%% 

hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
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
xlabel('EX-IN horizontal extent (µm)');
%stats against 0
[p,h] = signrank(diffL5fr_span);
%% Overview of centroid x any with respect to soma for EX and IN, 
a=iviv_cells;
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,0,[],[]);
%% 
% Creating correlation matrix OSI, DSI, ODI
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0 & abs(od_out_iviv(:,3))>0); 

% binodir=reshape([str(a).Dir],[2,length(a)])';
%  binodir_delta=binodir(:,1)-binodir(:,2);
% idx_bi=find(abs(od_out_iviv(a,3))<0.1)
% id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
% a(idx_bi(id_ov))=[];
% a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0 & abs(od_out_iviv(:,3))>0); 
% binodir=reshape([str(a).Ori],[2,length(a)])';
%  binodir_delta=binodir(:,1)-binodir(:,2);
% idx_bi=find(abs(od_out_iviv(a,3))<0.1)
% id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
% a(idx_bi(id_ov))=[];
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2) span(a,:) tot_inputL23(a,:) tot_inputL4(a,:) tot_inputL5(a,:) od_out_iviv(a,[1 2])]; 
[R1,P1]=corrcoef(com,'rows','pairwise');
G=correlation_matrix(com,0);title('');xticks([1:1:30]);yticks([1:1:30]);

%% Subsample sigma based on quality of the fit
a=[];
a=find(r_sq>0.3);
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2) span(a,:) tot_inputL23(a,:) tot_inputL4(a,:) tot_inputL5(a,:) od_out_iviv(a,[7])]; 
G2=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
[R2,P2]=corrcoef(com,'rows','pairwise');
%% 
a=[];
a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0 & abs(od_out_iviv(:,3))>0); 
% binodir=reshape([str(a).Ori],[2,length(a)])';
%  binodir_delta=binodir(:,1)-binodir(:,2);
% idx_bi=find(abs(od_out_iviv(a,3))<0.1)
% id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
% a(idx_bi(id_ov))=[];
com=[];com=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2) span(a,:) tot_inputL23(a,:) tot_inputL4(a,:) tot_inputL5(a,:) od_out_iviv(a,[6])]; 
[R3,P3]=corrcoef(com,'rows','pairwise');
G2=correlation_matrix(com,0);title('');xticks([1:1:16]);yticks([1:1:16]);
%% 

fG=[];fG2=[];
fG=[R1([19 20],1:18)];
fG2=[R2(19,1:18)];
fG3=[R3(19,1:18)];
fGf=[fG; fG2;fG3];
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 200])
imagesc(fGf);c=colorbar;pos = get(c,'Position');
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);xticks([1:1:18]);yticks([1:1:12])
yticklabels({'gOSI','gDSI','TW','Ca_{peak}'});set(gca,'FontSize',10);
xticklabels({'L2/3fr','L2/3fr','L4fr','L4fr','L5fr','L5fr','HEL2/3 ','HEL2/3','HEL4 ','HEL4','HEL5 ','HEL5','TL2/3 ','TL2/3','TL4 ','TL4','TL5 ','TL5'});xtickangle(45);set(gca,'FontSize',10);

%% 
%% Circular correlation for ORI
a=[];
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0 & abs(od_out_iviv(:,3))>0); 
binodir=reshape([str(a).Ori],[2,length(a)])';
binodir_delta=binodir(:,1)-binodir(:,2);
idx_bi=find(abs(od_out_iviv(a,3))<0.1)
id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
a(idx_bi(id_ov))=[];
rho1=[];
pval1=[];
par_c=[];
 par_c=[L23fr(a,1)  L23fr(a,2) L4fr(a,1)  L4fr(a,2) L5fr(a,1) L5fr(a,2) span(a,1) span(a,2) span(a,4)  tot_inputL23(a,1) tot_inputL23(a,2) tot_inputL4(a,1)...
     ]
for i=1:size(par_c,2)
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(a,4)), par_c(:,i))
end
%% 

tempy=[];
tempy=find(~isnan(span(a,3)));
[rhx phx] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),span(a(tempy),3))
tempy=[];
tempy=find(~isnan(span(a,5)));
[rhx2 phx2] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),span(a(tempy),5))
tempy=[];
tempy=find(~isnan(span(a,6)));
[rhx3 phx2] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),span(a(tempy),6))
tempy=[];
tempy=find(~isnan(tot_inputL4(a,2)));
[rhx4 phx4] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),tot_inputL4(a(tempy),2))
tempy=[];
tempy=find(~isnan(tot_inputL5(a,1)));
[rhx5 phx5] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),tot_inputL5(a(tempy),1))
tempy=[];
tempy=find(~isnan(tot_inputL5(a,2)));
[rhx6 phx6] = circ_corrcl(deg2rad(od_out_iviv(a(tempy),4)),tot_inputL5(a(tempy),2))
%% 
%% 

%% Absolute Cx relative to axis of displacement
% close all;
% %a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
% a=[];
% a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25 & abs(od_out_iviv(:,3))>0); 
% binodir=reshape([str(a).Dir],[2,length(a)])';
%  binodir_delta=binodir(:,1)-binodir(:,2);
% % 
% idx_bi=find(abs(od_out_iviv(a,3))<0.1)
% id_ov=find(abs(binodir_delta(find(abs(od_out_iviv(a,3))<0.1)))>90)
% a(idx_bi(id_ov))=[];
% fig3= figure;set(fig3, 'Name', 'Input distribution');set(fig3, 'Position', [200, 600, 200, 200]);set(gcf,'color','w');
% scatter(od_out_iviv(a,4),(ang_inL23(a,3)-ang_inL23(a,1))*69,12,'bo','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
% %hold on;scatter(od_out_iviv(a,4)-45,(ang_inL23(a,3)-ang_inL23(a,1))*69,12,'ro','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
% [rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), (ang_inL23(a,3)-ang_inL23(a,1))*69);
% text(-20,90,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','b')
% hold on
% scatter(od_out_iviv(a,4)-45,(ang_exL23(a,3)-ang_exL23(a,1))*69,12,'ro','filled');xlabel('\Delta  Angle from aligned');ylabel('|C_{x} (µm)|');
% [rho1 pval1] = circ_corrcl(deg2rad(od_out_iviv(a,4)), (ang_exL23(a,3)-ang_exL23(a,1))*69);
% text(-20,100,['cc=' num2str(round(rho1,3)) ' ' 'p<0.001'],'Color','r')
% set(gca,'FontSize',10);
% ylim([-10 100])
%% Morpholgy with relation to Orientation prference 
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
%% 
df=[morph_parameters(:,9) morph_parameters(:,10)];
db=[morph_parameters(:,19) morph_parameters(:,20)];
a=[];
a=morph_res_sub
 par_c=[morph_parameters(a,2) df(a,1) dis_s(a)'  morph_parameters(a,4)  max_s(a)' ...
    morph_parameters(a,12) db(a,2) dis_s_ba(a)'  morph_parameters(a,14) max_s_ba(a)']
%% Circular correlation for ORI
a=[];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
par_c=[];rho1=[];pval1=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),4)), par_c(:,i))
end
%% 

%%  barplot
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
s2=[g3'];
s1=[g4'];
par=df(a,1);
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL','NAL'});ylabel('Max horizontal extent apical');set(gca,'FontSize',10)
%title('Apical tree')
%% 
mcell=find([str(:).iviv]==1 & morph_cells==1);
morph_sel=[];
morph_sel=[morph_parameters(mcell,[2]) morph_parameters(mcell,[9]) dis_s(mcell)' morph_parameters(mcell,[4]) max_s(mcell)' ...
    morph_parameters(mcell,[12]) morph_parameters(mcell,[19]) dis_s_ba(mcell)' morph_parameters(mcell,[14]) max_s_ba(mcell)'];
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
%% 
%% Plot correlation matrix, PANEL F
df=[morph_parameters(morph_res_sub,9) morph_parameters(morph_res_sub,10)];
db=[morph_parameters(morph_res_sub,19) morph_parameters(morph_res_sub,20)];
com=[];com=[morph_parameters(morph_res_sub,2) nanmax(df,[],2) dis_s(morph_res_sub)'  morph_parameters(morph_res_sub,4)  max_s(morph_res_sub)' ...
    morph_parameters(morph_res_sub,12) nanmax(db,[],2) dis_s_ba(morph_res_sub)'  morph_parameters(morph_res_sub,14) max_s_ba(morph_res_sub)'  od_out_iviv(morph_res_sub,[1 2 3 6])]
G1=correlation_matrix(com,0);
%% subsample TW
df=[];db=[];
df=[morph_parameters(morph_res_sub2,9) morph_parameters(morph_res_sub2,10)];
db=[morph_parameters(morph_res_sub2,19) morph_parameters(morph_res_sub2,20)];
com=[];com=[morph_parameters(morph_res_sub2,2) nanmax(df,[],2) dis_s(morph_res_sub2)'  morph_parameters(morph_res_sub2,4)  max_s(morph_res_sub2)' ...
    morph_parameters(morph_res_sub2,12) nanmax(db,[],2) dis_s_ba(morph_res_sub2)'  morph_parameters(morph_res_sub2,14) max_s_ba(morph_res_sub2)'  od_out_iviv(morph_res_sub2,[7])]
G2=[];
G2=correlation_matrix(com,0);
%% Circular correlation for ORI
a=[];
a=find(od_out_iviv(morph_res_sub,1)>0.25 & r_sq(morph_res_sub)>0.3) ; 
par_c=[];rho1=[];pval1=[];
 par_c=[morph_parameters(morph_res_sub(a),2) nanmax(df(a),[],2) dis_s(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),4)  max_s(morph_res_sub(a))' ...
    morph_parameters(morph_res_sub(a),12) nanmax(db(a),[],2) dis_s_ba(morph_res_sub(a))'  morph_parameters(morph_res_sub(a),14) max_s_ba(morph_res_sub(a))']
for i=1:10
    [rho1(i) pval1(i)] = circ_corrcl(deg2rad(od_out_iviv(morph_res_sub(a),4)), par_c(:,i))
end
%% 

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
%% only signifcant from circ
c_cor=[rho1;rho2]
c_pva=[pval1;pval2]
c_cor_e=c_cor;
  m=c_pva<0.05;
c_cor_e(m==0)=m(m==0);
%% 
Gf=G1(11:14,1:10)
tG=[Gf;G2(11:end,1:10);c_cor_e;G1(14,1:10)]
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [400, 600, 500, 325]);imagesc(tG(1:7,:));c=colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);caxis([-1 1]);
xticks([1:1:12]);yticks([1:1:7]) 
xticklabels({'Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing','Total Length','Max extent','Dis peak branch','Nr. branch points','Peak number crossing'});xtickangle(45);set(gca,'FontSize',12)
yticklabels({'gOSI','gDSI','ODI','Ca_{peak}','TW','ORI','DIR','Pial depth'});set(gca,'FontSize',12)
c.Label.String = 'r';set(gca,'FontSize',10); c.Ticks=[-1:0.5:1]; set(gca,'FontSize',10);

%% 
[max_s_sub dis_s_sub max_s_ba_sub dis_s_ba_sub]=sholl_analysis(zz,mcell,1);
%% 
set(gcf,'color','w');
box off;
legend('Apical','Basal');legend boxoff;ylabel('Nr. dendritic crossings');xlabel('Distance from Soma (µm)');set(gca,'FontSize',10);
%% 

%% Read out orientation for barplots
a=[];
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0); 
sector=50;
midpoint=135;
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
s2=[g3'];
s1=[g4'];
%par=span(a,1)*69;
par=span(a,4)*69;
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);hold on;xticklabels({'AL','NAL'});ylabel('L2/3 Horizontal extent (µm)');set(gca,'FontSize',10)