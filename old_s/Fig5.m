%Figure 5
%% 
%Overview of centroid x and y with respect to soma for EX and IN and OSI colurcoded, Panel B
a=find(od_out_iviv(:,1)>0.25);
fe=od_out_iviv(a,4);
centroid_plot(a,ang_exL23,ang_exL4,ang_exL5,ang_inL23,ang_inL4,ang_inL5,1,fe,{'ORI'});
%% 
%% Rolling average  Cy
parameter_vector = abs(ang_exL23(:,4)*69-ang_exL23(:,2)*69)
parameter_vector2 = abs(ang_inL23(:,4)*69-ang_inL23(:,2)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L2/3 (�m)','Color','k');
  ylim([5 70]);yticks([10:30:70]);
 xlim([0 180]);xticks([0:10:180]);
set(gca,'FontSize',10);
%% Circular correlation for Cy EX and IN
a=find(od_out_iviv(:,1)>0.25);  
[rho pval] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs(ang_exL23(a,4)*69-ang_exL23(a,2)*69))
a=find(od_out_iviv(:,1)>0.25);  
%% Barplots Cy ex and EX
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_exL23(a,4)*69-ang_exL23(a,2)*69)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (�m)','Color','k'); ;set(gca,'FontSize',10); 
%% Barplots Cy ex and IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL23(a,4)*69-ang_inL23(a,2)*69)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<170);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (�m)','Color','k'); ;set(gca,'FontSize',10);
%% Pial depth
parameter_vector = pia_input
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,0)
ylabel('Pial depth (�m)','Color','k');
ylim([170 270]);yticks([170:50:280]);
xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
set(gca,'Ydir','reverse')
%% Pial depth barplot
a=find(od_out_iviv(:,1)>0.25);  
par=pia_input(a)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<170);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Pial depth (�m)','Color','k'); ;set(gca,'FontSize',10);

%% 
%% Rolling average  Cy
parameter_vector = ang_exL23(:,4)*69
parameter_vector2 = ang_inL23(:,4)*69
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L2/3 (�m)','Color','k');
ylim([170 270]);yticks([170:50:280]);
  xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
set(gca,'Ydir','reverse')
%% Cy raw ex
a=find(od_out_iviv(:,1)>0.25);  
par=ang_exL23(a,4)*69
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<170);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy raw (�m)','Color','k'); ;set(gca,'FontSize',10);
%% Cy raw ex
a=find(od_out_iviv(:,1)>0.25);  
par=ang_inL23(a,4)*69
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<170);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy raw (�m)','Color','k'); ;set(gca,'FontSize',10);
%% 
%% L4 fraction
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,1)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','k'); ;set(gca,'FontSize',10); 
ylim([0 0.6]);yticks([0:0.3:0.6]);
%% %% L4 fraction
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,2)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','k'); ;set(gca,'FontSize',10); 
ylim([0 0.6]);yticks([0:0.3:0.6]);
%% 
a=find(od_out_iviv(:,1)>0.25);
 corr_plot((ang_inL23(a,4)*69-ang_inL23(a,2)*69)*-1,pia_input(a),od_out_iviv(a,4),{'','',''});
% xlabel('Cy L23');ylabel('Pial depth (�m)','Color','b');set(gca,'Ydir','reverse')
%corr_plot(ang_inL23(a,4)*69,pia_input(a),od_out_iviv(a,4),{'','',''});
xlabel('Cy L23');ylabel('Pial depth (�m)','Color','b');set(gca,'Ydir','reverse')
%% 
a=find(od_out_iviv(:,1)>0.25);
 corr_plot(ang_inL23(a,4),pia_input(a),od_out_iviv(a,4),{'','',''});
% xlabel('Cy L23');ylabel('Pial depth (�m)','Color','b');set(gca,'Ydir','reverse')
%corr_plot(ang_inL23(a,4)*69,pia_input(a),od_out_iviv(a,4),{'','',''});
xlabel('Cy L23');ylabel('Pial depth (�m)','Color','b');set(gca,'Ydir','reverse')

%% 
figure;set(gcf,'color','w')
q=quiver(ang_inL23(a,1)*69,ang_inL23(a,2)*69,ones(length(a),1)*69*2,(ang_inL23(a,4)-ang_inL23(a,2))*69,0);
 q.Color='black';q.MaxHeadSize=0.1;q.LineWidth=0.1;
hold on;scatter(ang_inL23(a,1)*69,ang_inL23(a,2)*69,60,od_out_iviv(a,4),'filled','^');cmap=phasemap;colormap(cmap)
set(gca,'Ydir','reverse');box off;ylabel('Pial depth (�m)');
%% 
%% 
parameter_vector = abs(ang_exL4(:,4)*69-ang_exL4(:,2)*69)
parameter_vector2 = abs(ang_inL4(:,4)*69-ang_inL4(:,2)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L4 (�m)','Color','k');
%   ylim([5 70]);yticks([10:30:70]);
%  xlim([0 180]);xticks([0:10:180]);
set(gca,'FontSize',10);
%% 

a=find(od_out_iviv(:,1)>0.25);  
par=ang_exL4(a,4)*69
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (�m)','Color','k'); ;set(gca,'FontSize',10); 

%% 
figure
a=find(od_out_iviv(:,1)>0.25);  
h4 = histogram((ang_exL23(a,4)*69-ang_exL23(a,2)*69)-(ang_inL23(a,4)*69-ang_inL23(a,2)*69),8);h4.EdgeColor = 'k';h4.FaceColor = 'm';hold on;
box off;xlabel('Cy delta EX-IN')
%% 
%[aa bb]=intersect(find(~isnan(ang_exL23(:,3)*69-ang_exL23(:,1)*69)),find(~isnan(ang_inL23(:,3)*69-ang_inL23(:,1)*69)));
[p1]=signrank(abs(ang_exL23(a,4)*69-ang_exL23(a,2)*69),abs(ang_inL23(a,4)*69-ang_inL23(a,2)*69))
%% Cx part

%% Rolling average  Cx
parameter_vector = abs(ang_exL23(:,3)*69-ang_exL23(:,1)*69)
parameter_vector2 = abs(ang_inL23(:,3)*69-ang_inL23(:,1)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cx L2/3 (�m)','Color','k');
  ylim([5 70]);yticks([10:30:70]);
 xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
%% 
%% Barplots Cx ex and EX
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);  
par=abs(ang_exL23(a,3)*69-ang_exL23(a,1)*69)
g1=find(od_out_iviv(a,4)>0 & od_out_iviv(a,4)<50) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylim([0 110]);yticks([0:55:110]);
ylabel('Cx L2/3 (�m)','Color','k'); ;set(gca,'FontSize',10); 
%% Barplots Cx ex and IN
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);  
par=abs(ang_inL23(a,3)*69-ang_inL23(a,1)*69)
g1=find(od_out_iviv(a,4)>0 & od_out_iviv(a,4)<50) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylim([0 110]);yticks([0:55:110])
ylabel('Cx L2/3 (�m)','Color','k') ;set(gca,'FontSize',10);
%% No corrlaion between Cx and pia
a=find(od_out_iviv(:,1)>0.25);
figure;set(gcf,'color','w')
 scatter(ang_exL23(a,3)*69-ang_exL23(a,1)*69,pia_input(a),14,'r','filled');
 hold on;scatter(ang_inL23(a,3)*69-ang_inL23(a,1)*69,pia_input(a),14,'b','filled');
set(gca,'Ydir','reverse');box off;ylabel('Pial depth (�m)');
ylim([170 270]);yticks([170:50:280]);
xlim([-120 120]);xticks([-120:60:120]);xlabel('Cx L2/3 (�m)');set(gca,'FontSize',10); 
%% 


%% Non visual responsive
nvr=find([str(:).sftf_resp]==0 & [str(:).resp]==0);
nvra=find([str(:).resp]==0);
nvrb=find([str(:).sftf_resp]==1 & [str(:).resp]==0);
%nvr=find([str(:).resp]==0);
vr=find([str(:).resp]==1);
w1=find(od_out_iviv(:,4)>0 & od_out_iviv(:,4)<50 & od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);
w1a=find(od_out_iviv(:,5)>195 & od_out_iviv(:,5)<245 & od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25);
w1b=find(od_out_iviv(:,5)>15 & od_out_iviv(:,5)<65 & od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.25);
w2=find(od_out_iviv(:,4)>100 & od_out_iviv(:,4)<150 & od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0);
%%  
par=abs(ang_inL23(:,3)*69-ang_inL23(:,1)*69)
g1=nvra
g2=w1
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'non resp.',''});xtickangle(45);
%ylim([0 110]);yticks([0:55:110])
ylabel('Cx L2/3 (�m)','Color','k') ;set(gca,'FontSize',10);
%% 
par1=abs(ang_inL23(:,3)*69-ang_inL23(:,1)*69)
%par1=span(:,6);
%par1=ang_inL5(:,3)*69-ang_inL5(:,1)*69
%par1=ang_inL5(:,8)*69
data= vertcat(par1(nvra),par1(w1),par1(w2));
groups_idx=vertcat(ones(length(nvra),1)*1,ones(length(w1),1)*2,ones(length(w2),1)*3)
[statsout] = barplot_sw(data,groups_idx,{'','Cd length'});

 %% 
 g1=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.15 & od_out_iviv(:,5)>270 & od_out_iviv(:,5)<330); 
  g2=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.15 & od_out_iviv(:,5)>90 & od_out_iviv(:,5)<150); 
 figure;set(gcf,'color','w')
 scatter(ang_exL23(g1,3)-ang_inL23(g1,3),pia_input(g1),'r','filled');
 hold on;
 scatter(ang_exL23(g2,3)-ang_inL23(g2,3),pia_input(g2),'b','filled');
%% 
g1=find(od_out_iviv(:,1)>0 & od_out_iviv(:,2)>0.15 & od_out_iviv(:,5)>270 & od_out_iviv(:,5)<330); 
figure;histogram(od_out_iviv(g1,5))
figure;scatter(ang_inL23(g1,3)*69-ang_inL23(g1,1)*69,ones(length(g1),1),20,od_out_iviv(g1,5),'filled')
%% 

%% 


[rho pval] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs(ang_inL23(a,4)*69-ang_inL23(a,2)*69))
a=find(od_out_iviv(:,1)>0.25); 
[rho pval] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs(ang_exL23(a,3)*69-ang_exL23(a,1)*69))
a=find(od_out_iviv(:,1)>0.25);  
[rho pval] = circ_corrcl(deg2rad(od_out_iviv(a,4)), abs(ang_inL23(a,3)*69-ang_inL23(a,1)*69))