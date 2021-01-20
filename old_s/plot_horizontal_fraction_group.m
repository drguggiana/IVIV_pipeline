function statsout=plot_horizontal_fraction_group(str,od_out_iviv,L23h,L4h,L5h)
%% Plot horizontal fraction per layer per direction group
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
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color','r');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color','b');
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
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s1,:),1)),'--k');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'--k');
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');
xlim([1 16]);ylim([-0.35 0.35]);yticks(-0.35:0.35:0.35);
xticks(1:5:16);xticklabels({'-552','-138','138','552'})

%text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
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
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color','r');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color','b');
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
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s1,:),1)),'--k');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'--k');
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
xlim([1 16])
%hold on;text(2.5,0.55,'L4','FontSize',12);
xlabel('Horizontal position (µm)')
xticks(1:5:16);xticklabels({'-552','-138','138','552'})
ylim([-0.35 0.35]);yticks(-0.35:0.35:0.35);
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
mexp=errorbar(nanmean( flip_hl),nanstd( flip_hl)/sqrt(size( flip_hl,1)),'Color','r');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean( flip_hl_in)*-1,nanstd( flip_hl_in)/sqrt(size( flip_hl_in,1))*-1,'Color','b');
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
mexp=errorbar(nanmean(frac_h(s2,1:16)),nanstd(frac_h(s2,1:16))/sqrt(size(frac_h(s2,:),1)),'--k');
hold on;mexp.CapSize=3;
hold on;
mexp=errorbar(nanmean(frac_h(s2,17:end))*-1,nanstd(frac_h(s2,17:end))/sqrt(size(frac_h(s2,:),1))*-1,'--k');
hold on;mexp.CapSize=3;
hold on;line([8.5 8.5], [-0.4 0.4],'Color','k','LineStyle','--');text(2.5,0.4,'M');text(13,0.4,'L');set(gca,'FontSize',10);
%hold on;text(2.5,0.55,'L5','FontSize',12);
xlim([1 16]);
xticks(1:5:16);xticklabels({'-552','-138','138','552'})
ylim([-0.35 0.35]);yticks(-0.35:0.35:0.35);

%% 
%% Statistics
for i=1:16;
par=[];
par=L23h(a,i);
t=0;
if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par(s1),par(s2));
    stats1(i)=p;
elseif t==1
    [p k]=ranksum(par(s1),par(s2),'tail','right');
    stats1(i)=p;
else t==2
    [p k]=ranksum(par(s1),par(s2),'tail','left');
    stats1(i)=p;
end
end

for i=1:16;
par=[];
par=L4h(a,i);
t=0;
if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par(s1),par(s2));
    stats2(i)=p;
elseif t==1
    [p k]=ranksum(par(s1),par(s2),'tail','right');
    stats2(i)=p;
else t==2
    [p k]=ranksum(par(s1),par(s2),'tail','left');
    stats2(i)=p;
end
end

for i=1:16;
par=[];
par=L5h(a,i);
t=0;
if t==0
%    [k p]=ttest2(par(g1),par(g2));
%    statsout=p;
[p k]=ranksum(par(s1),par(s2));
    stats3(i)=p;
elseif t==1
    [p k]=ranksum(par(s1),par(s2),'tail','right');
    stats3(i)=p;
else t==2
    [p k]=ranksum(par(s1),par(s2),'tail','left');
    stats3(i)=p;
end
end

statsout=[stats1; stats2; stats3]';
end