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
ylabel('Cy L2/3 (µm)','Color','k');
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
ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10); 
%% Barplots Cy ex and IN
a=find(od_out_iviv(:,1)>0.25);  
par=abs(ang_inL23(a,4)*69-ang_inL23(a,2)*69)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<170);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10);
%% Pial depth
parameter_vector = pia_input
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,0)
ylabel('Pial depth (µm)','Color','k');
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
ylabel('Pial depth (µm)','Color','k'); ;set(gca,'FontSize',10);

%% 
%% Rolling average  Cy
parameter_vector = ang_exL23(:,4)*69
parameter_vector2 = ang_inL23(:,4)*69
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L2/3 (µm)','Color','k');
ylim([170 270]);yticks([170:50:280]);
  xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
set(gca,'Ydir','reverse')
%% 
%% L4 fraction
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,1)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','r'); ;set(gca,'FontSize',10); 
%% %% L4 fraction
a=find(od_out_iviv(:,1)>0.25);  
par=L4fr(a,2)
g1=find(od_out_iviv(a,4)>20 & od_out_iviv(a,4)<90) ;
g2=find(od_out_iviv(a,4)>115 & od_out_iviv(a,4)<185);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('L4fr','Color','b'); ;set(gca,'FontSize',10); 
%% 
a=find(od_out_iviv(:,1)>0.25);
 corr_plot((ang_inL23(a,4)*69-ang_inL23(a,2)*69)*-1,pia_input(a),od_out_iviv(a,4),{'','',''});
% xlabel('Cy L23');ylabel('Pial depth (µm)','Color','b');set(gca,'Ydir','reverse')
%corr_plot(ang_inL23(a,4)*69,pia_input(a),od_out_iviv(a,4),{'','',''});
xlabel('Cy L23');ylabel('Pial depth (µm)','Color','b');set(gca,'Ydir','reverse')
%% 
a=find(od_out_iviv(:,1)>0.25);
 corr_plot(ang_inL23(a,4),pia_input(a),od_out_iviv(a,4),{'','',''});
% xlabel('Cy L23');ylabel('Pial depth (µm)','Color','b');set(gca,'Ydir','reverse')
%corr_plot(ang_inL23(a,4)*69,pia_input(a),od_out_iviv(a,4),{'','',''});
xlabel('Cy L23');ylabel('Pial depth (µm)','Color','b');set(gca,'Ydir','reverse')

%% 
figure;set(gcf,'color','w')
q=quiver(ang_inL23(a,1)*69,ang_inL23(a,2)*69,ones(length(a),1)*69*2,(ang_inL23(a,4)-ang_inL23(a,2))*69,0);
 q.Color='black';q.MaxHeadSize=0.1;q.LineWidth=0.1;
hold on;scatter(ang_inL23(a,1)*69,ang_inL23(a,2)*69,60,od_out_iviv(a,4),'filled','^');cmap=phasemap;colormap(cmap)
set(gca,'Ydir','reverse');box off;ylabel('Pial depth (µm)');
%% 
%% 
parameter_vector = abs(ang_exL4(:,4)*69-ang_exL4(:,2)*69)
parameter_vector2 = abs(ang_inL4(:,4)*69-ang_inL4(:,2)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cy L4 (µm)','Color','k');
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
ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10); 
%% Cx part

%% Rolling average  Cx
parameter_vector = abs(ang_exL23(:,3)*69-ang_exL23(:,1)*69)
parameter_vector2 = abs(ang_inL23(:,3)*69-ang_inL23(:,1)*69)
rolling_avg_display(str,parameter_vector,parameter_vector2,45,1,1)
ylabel('Cx L2/3 (µm)','Color','k');
  ylim([5 70]);yticks([10:30:70]);
 xlim([0 180]);xticks([0:45:180]);
set(gca,'FontSize',10);
%% 
%% Barplots Cx ex and EX
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0.15);  
par=abs(ang_exL23(a,3)*69-ang_exL23(a,1)*69)
g1=find(od_out_iviv(a,4)>0 & od_out_iviv(a,4)<50) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (µm)','Color','k'); ;set(gca,'FontSize',10); 
%% Barplots Cx ex and IN
a=find(od_out_iviv(:,1)>0.25 & od_out_iviv(:,2)>0.15);  
par=abs(ang_inL23(a,3)*69-ang_inL23(a,1)*69)
g1=find(od_out_iviv(a,4)>0 & od_out_iviv(a,4)<50) ;
g2=find(od_out_iviv(a,4)>100 & od_out_iviv(a,4)<150);
[statsout]=dual_barplot(par,g1,g2,1);xticks([1:1:2]);
xticklabels({'',''});xtickangle(45);
ylabel('Cy L2/3 (µm)','Color','k') ;set(gca,'FontSize',10);
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