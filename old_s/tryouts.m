
a=cell_idx
pv = (ang_exL23(a,3)-ang_exL23(a,1))*69
corr_plot(od_out_iviv(a,5),score_in(:,2),od_out_iviv(a,4),{'Orientation along slice','IN',''});



%% 
a=cell_idx
tt=ones(1,147)*NaN;
tt(a)=score_in(:,2);
tt2(a)=score_ex(:,3);

parameter_vector = tt
parameter_vector2 = tt2
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
%% 

s1a=10;s1b=65;s2a=100;s2b=155;
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0.25);  
par=abs(tt(a))
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);

%% 

mm=od_out_iviv(a,5)
par1=[tt(a)]
par2=[tt2(a)]

rt=[par1(g1)' par2(g1)']

t1=1:6
t2=7:12
[statsout]=dual_barplot(rt,t1,t2,0);xticks([1:1:2]);

%% 
paired_plot(rt,1,{'b','r'})
%% 
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0.2);  
pv = (ang_exL23(a,3)-ang_exL23(a,1))*69
corr_plot(od_out_iviv(a,4),abs(pv),od_out_iviv(a,2),{'Orientation along slice','IN',''});
%% 

s1a=290;s1b=340;s2a=110;s2b=165;
a=find(od_out_iviv(:,1)>0.2 & od_out_iviv(:,2)>0.2);  
par=abs((ang_inL23(a,3)-ang_inL23(a,1))*69)
g1=find(od_out_iviv(a,5)>s1a & od_out_iviv(a,5)<s1b) ;
g2=find(od_out_iviv(a,4)>s2a & od_out_iviv(a,4)<s2b);
[statsout]=dual_barplot(par,g1,g2,0);xticks([1:1:2]);
%% 

mm=od_out_iviv(a,5)
par1=[tt(a)]
par2=[tt2(a)]

rt=[par1(g1)' par2(g1)']


paired_plot(rt,1,{'b','r'})
%% 
s1a=10;s1b=60;s2a=110;s2b=165;

a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)<0.25);  
par=(ang_inL23(a,3)-ang_inL23(a,1))*69
g1=find(od_out_iviv(a,4)>s1a & od_out_iviv(a,4)<s1b) ;
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

s1=[g1' g2'];
s2=[g3' g4'];
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
par(g1)=par(g1)*-1
[statsout]=dual_barplot(par,s1,s2,1);xticks([1:1:2]);
%% 
tempr=od_out_iviv(a,3)
figure;set(gcf,'color','w');scatter(par(s1),abs(tempr(s1)),'m','filled');
ylabel('|ODI|');xlabel('Centroid offset')
%hold on;scatter(par(s2),abs(tempr(s2)),'k','filled');
%% 




 

corr_plot((ang_exL23(a,3)-ang_exL23(a,1))*69,(ang_inL23(a,3)-ang_inL23(a,1))*69,[],{'Orientation along slice','IN',''});
%% 

mm=od_out_iviv(a,5)
par=(ang_exL23(a,3)-ang_exL23(a,1))*69;
[statsout]=dual_barplot(par,g1,g2,2);xticks([1:1:2]);

%% 
idxb = discretize(mm,4);
for i=1:max(idxb)
avg_bin(:,i)=nanmean(par(idxb==i));
avg_sem(:,i)=nanstd(par(idxb==i))/sqrt(length(par(idxb==i)));
end


figure;errorbar(avg_bin,avg_sem)

%% 



parameter_vector = (ang_exL23(:,3)-ang_exL23(:,1))*69;
parameter_vector2 = (ang_inL23(:,3)-ang_inL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1)
%% 
ylabel('Centroid offset')

%% 
%% Calculate difference between ex and in (ex-in) medial and lateral for whole map
%mean
frh_medial=[nanmean(frh(:,1:8),2) nanmean(frh(:,17:24),2)];
frh_lateral=[nanmean(frh(:,9:16),2) nanmean(frh(:,25:end),2)];

% %sum
frh_medial_s=[sum(frh(:,1:8),2) sum(frh(:,17:24),2)];
frh_lateral_s=[sum(frh(:,9:16),2) sum(frh(:,25:end),2)];
%% 

parameter_vector = frh_lateral_s(:,1)-0.5
parameter_vector2 = frh_lateral_s(:,2)-0.5
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1);

%% 
ylabel('Delta Fraction lateral')
%% 

par=frh_medial_s(a,1)-0.5
[statsout]=dual_barplot(abs(par),s1,s2,1);xticks([1:1:2]);
%% 
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
%% 

corr_plot(frh_medial_s(a,2),(ang_exL23(a,3)-ang_exL23(a,1))*69,od_out_iviv(a,2),{'Orientation along slice','IN',''});

%% Span in rolling average
parameter_vector = span(:,1)*69;
parameter_vector2 = span(:,4)*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1);
ylabel('Horizontal extent')
%% Span barplot
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

s1=[g3'];
s2=[g4'];
par=span(a,4)*69;
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);
%% Max pixel
mx_all=[];
mx_all_e=[];
for i=1:147
my=[];mx=[];
temp=[];
temp=in_map(3:5,1:16,i);
[my mx]=find(temp==max(max(temp)));
mx_all(i)=mx(1);
end
for i=1:147
my=[];mx=[];
temp=[];
temp=ex_map(3:5,1:16,i);
[my mx]=find(temp==max(max(temp)));
mx_all_e(i)=mx(1);
end
%% Max in rolling average
parameter_vector = (mx_all_e'-ang_exL23(:,1))*69;
parameter_vector2 = (mx_all'-ang_exL23(:,1))*69;
rolling_avg_display(str,parameter_vector,parameter_vector2,65,0,1)

%% 

ylabel('Max pixel offset')
%% 
a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
corr_plot((mx_all(a)'-ang_inL23(a,1))*69,(ang_inL23(a,3)-ang_exL23(a,1))*69,od_out_iviv(a,5),{'Orientation along slice','IN',''});

figure;histogram((mx_all(a)'-ang_inL23(a,1))*69-(ang_exL23(a,3)-ang_exL23(a,1))*69)

%% 

[stats_g] = display_inputs([frv],[frh],frv(:,1:16)-frv(:,17:end),frh(:,1:16)-frh(:,17:end),[]);
%% 

a=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25); 
sector=60;
midpoint=135;

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
%% 

fr_s=[];
fr_vs=[];
fr_s=frh(a,:);
fr_vs=frv(a,:);

fr_s(g1,:)=flip(fr_s(g1,:));
fr_vs(g1,:)=flip(fr_vs(g1,:));
%% 


%% 
[stats_g] = display_inputs([fr_vs(s1,:)],[fr_s(s1,:)],fr_vs(s1,1:16)-fr_vs(s1,17:end),fr_s(s1,1:16)-fr_s(s1,17:end),[]);

%% 

[stats_g] = display_inputs([fr_vs(s2,:)],[fr_s(s2,:)],fr_vs(s2,1:16)-fr_vs(s2,17:end),fr_s(s2,1:16)-fr_s(s2,17:end),[]);

%% 
su_d=[fr_s(s1,:); fr_s(s2,:)]
su_vd=[fr_vs(s1,:); fr_vs(s2,:)]

group=[ones(10,1); ones(18,1)*2]
[stats_g] = display_inputs([su_vd(:,:)],[su_d(:,:)],su_vd(:,1:16)-su_vd(:,17:end),su_d(:,1:16)-su_d(:,17:end),group');
%% 
frac_v=[];
frac_h=[];
frac_v=fr_vs(s1,:);
frac_h=fr_s(s1,:);
figure;
for i=1:size(frac_v,1)
exp=plot(frac_h(i,1:16)','--r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
exp=plot(frac_h(i,17:end)'*-1,'--b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
hold on
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-1 1]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,1,'Medial');text(10,1,'Lateral');set(gca,'FontSize',10);


hold on;

frac_v=[];
frac_h=[];
frac_v=fr_vs(s2,:);
frac_h=fr_s(s2,:);
for i=1:size(frac_v,1)
exp=plot(frac_h(i,1:16)','-r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
exp=plot(frac_h(i,17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'m');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'m');
mexp.CapSize=3;ylim([-0.35 0.35]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,1,'Medial');text(10,1,'Lateral');set(gca,'FontSize',10);

%% 
frac_diffh=[];
frac_diffh=fr_s(s1,1:16)-fr_s(s1,17:end)

figure;
mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'k');
mexp.CapSize=3;
ylim([-0.05 0.05]);
yticks([-0.05:0.05:0.05]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
hold on;text(5,0.105,'\Delta EX - IN','Color','k');
hold on;line([1 16], [0 0],'Color','k','LineStyle','--')

hold on
frac_diffh=[];
frac_diffh=fr_s(s2,1:16)-fr_s(s2,17:end)
mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'r');
mexp.CapSize=3;
ylim([-0.05 0.05]);
yticks([-0.05:0.05:0.05]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
hold on;text(5,0.105,'\Delta EX - IN','Color','k');
hold on;line([1 16], [0 0],'Color','k','LineStyle','--')

%% 

%% fractipn per layer
for i=1:147
my=[];mx=[];
temp=[];
temp=in_map(6:7,1:16,i);
mx_all(i,:)=sum(temp)/2;
temp=[];
end

for i=1:147
my=[];mx=[];
temp=[];
temp=ex_map(6:7,1:16,i);
mx_all_e(i,:)=sum(temp)/2;
temp=[];
end
%% 

fr_s=[];
fr_vs=[];

fr_s=[mx_all_e(a,:) mx_all(a,:)];

fr_s(g1,:)=flip(fr_s(g1,:));
%% 

frac_v=[];
frac_h=[];

frac_h=fr_s(s1,:);
figure;
set(gcf,'color','w');
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','--r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
exp=plot(frac_h(i,17:end)'*-1,'--b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
hold on
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'m');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'m');
mexp.CapSize=3;ylim([-3 3]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');text(2.5,1,'Medial');text(10,1,'Lateral');set(gca,'FontSize',10);


hold on;

frac_v=[];
frac_h=[];

frac_h=fr_s(s2,:);
for i=1:size(frac_h,1)
exp=plot(frac_h(i,1:16)','-r');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
exp=plot(frac_h(i,17:end)'*-1,'-b');ylabel('Fraction of total input');box off;
exp.Color(4) = 0.1;
end
mexp=errorbar(nanmean(frac_h(:,1:16)),nanstd(frac_h(:,1:16))/sqrt(size(frac_h,1)),'k');
hold on;mexp.CapSize=3;
mexp=errorbar(nanmean(frac_h(:,17:end))*-1,nanstd(frac_h(:,17:end))/sqrt(size(frac_h,1)),'k');
mexp.CapSize=3;ylim([-1 1]);yticks([-1:0.5:1]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-3 3],'Color','k','LineStyle','--');text(2.5,1,'Medial');text(10,1,'Lateral');set(gca,'FontSize',10);
%% 


frac_diffh=[];
frac_diffh=fr_s(s1,1:16)-fr_s(s1,17:end)

figure;set(gcf,'color','w');
mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'m');
mexp.CapSize=3;
ylim([-0.05 0.05]);
yticks([-0.05:0.05:0.05]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
hold on;text(5,0.105,'\Delta EX - IN','Color','k');
hold on;line([1 16], [0 0],'Color','k','LineStyle','--')

hold on
frac_diffh=[];
frac_diffh=fr_s(s2,1:16)-fr_s(s2,17:end)
mexp=errorbar(nanmean(frac_diffh(:,1:16)),nanstd(frac_diffh(:,1:16))/sqrt(size(frac_diffh,1)),'k');
mexp.CapSize=3;
ylim([-0.1 0.5]);
yticks([-0.05:0.05:0.05]);hold on;xlim([1 16]);xticks([1:5:16]);xticklabels({'1','6','11','16'});
hold on;line([8.5 8.5], [-1 1],'Color','k','LineStyle','--');set(gca,'FontSize',10);
hold on;text(5,0.105,'\Delta EX - IN','Color','k');
hold on;line([1 16], [0 0],'Color','k','LineStyle','--')
%% 



a=find(od_out_iviv(:,1)<1 & od_out_iviv(:,2)>0.25); 
b=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)<0.1); 



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
s2=b;
par1=(ang_inL23(a,3)-ang_inL23(a,1))*69;
par1(g1)=par1(g1)*-1
par2=(ang_inL23(b,3)-ang_inL23(b,1))*69;

par1(s1)

%% 

color_id={[0.7 0 0.4],[0 0.5 0.5]};
figure;set(gcf,'color','w');
for i=1:2
gr_m=[nanmean(par1(s1)) nanmean(par2)];  
gr_sem=[nanstd(par1(s1))/sqrt(length(par1(s1))) nanstd(par2)/sqrt(length(par2))];
hold on;
b=bar(i,gr_m(i));b.FaceColor=color_id{i};
hold on;
plot(ones(1,length(par1(s1))),par1(s1),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
hold on;
plot(ones(1,length(par2))*2,par2,'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
end


t=1
if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par1(s1),par2);
    statsout=p
elseif t==1
    [p k]=ranksum(par1(s1),par2,'tail','right');
    statsout=p
else t==2
    [p k]=ranksum(par(g1),par(g2),'tail','left');
    statsout=p
end

for i=1:2
hold on;
er=errorbar(i,gr_m(i),gr_sem(i));er.Color = [0 0 0];er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
end

ylabel('Centroid offset')
%% 

