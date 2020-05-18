
g1=find(od_out_iviv(:,1)>0.25 & ang_inL4(:,3)*69-ang_inL4(:,1)*69<150 & od_out_iviv(:,4)>10 & od_out_iviv(:,4)<60);
g2=find(od_out_iviv(:,1)>0.25 & ang_inL4(:,3)*69-ang_inL4(:,1)*69<150 & od_out_iviv(:,4)>100 & od_out_iviv(:,4)<150);
%% EX g1

par1=ang_exL23(g1,3)*69-ang_exL23(g1,1)*69
par2=ang_exL4(g1,3)*69-ang_exL4(g1,1)*69

[R1,P1,RLO1,RUP1]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)

par1=ang_exL23(g1,3)*69-ang_exL23(g1,1)*69
par2=ang_exL5(g1,3)*69-ang_exL5(g1,1)*69

[z2 p2]=corrcoef(par1,par2,'Rows','pairwise')

par1=ang_exL4(g1,3)*69-ang_exL4(g1,1)*69
par2=ang_exL5(g1,3)*69-ang_exL5(g1,1)*69

[z3 p3]=corrcoef(par1,par2,'Rows','pairwise')

%% IN g1

par1=ang_inL23(g1,3)*69-ang_inL23(g1,1)*69
par2=ang_inL4(g1,3)*69-ang_inL4(g1,1)*69

[R4,P4,RLO4,RUP4]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)

par1=ang_inL23(g1,3)*69-ang_inL23(g1,1)*69
par2=ang_inL5(g1,3)*69-ang_inL5(g1,1)*69

[z5 p5]=corrcoef(par1,par2,'Rows','pairwise')

par1=ang_inL4(g1,3)*69-ang_inL4(g1,1)*69
par2=ang_inL5(g1,3)*69-ang_inL5(g1,1)*69

[z6 p6]=corrcoef(par1,par2,'Rows','pairwise')

%% EX g2

par1=ang_exL23(g2,3)*69-ang_exL23(g2,1)*69
par2=ang_exL4(g2,3)*69-ang_exL4(g2,1)*69

[R7,P7,RLO7,RUP7]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)

par1=ang_exL23(g2,3)*69-ang_exL23(g2,1)*69
par2=ang_exL5(g2,3)*69-ang_exL5(g2,1)*69

[z8 p8]=corrcoef(par1,par2,'Rows','pairwise')

par1=ang_exL4(g2,3)*69-ang_exL4(g2,1)*69
par2=ang_exL5(g2,3)*69-ang_exL5(g2,1)*69

[z9 p9]=corrcoef(par1,par2,'Rows','pairwise')

%% IN g2
par1=ang_inL23(g2,3)*69-ang_inL23(g2,1)*69
par2=ang_inL4(g2,3)*69-ang_inL4(g2,1)*69

[R10,P10,RLO10,RUP10]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)

par1=ang_inL23(g2,3)*69-ang_inL23(g2,1)*69
par2=ang_inL5(g2,3)*69-ang_inL5(g2,1)*69

[z11 p11]=corrcoef(par1,par2,'Rows','pairwise')

par1=ang_inL4(g2,3)*69-ang_inL4(g2,1)*69
par2=ang_inL5(g2,3)*69-ang_inL5(g2,1)*69
%% 

par1=ang_exL23(nvra,3)*69-ang_exL23(nvra,1)*69
par2=ang_exL4(nvra,3)*69-ang_exL4(nvra,1)*69
[R13,P13,RLO13,RUP13]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)

par1=ang_inL23(nvra,3)*69-ang_inL23(nvra,1)*69
par2=ang_inL4(nvra,3)*69-ang_inL4(nvra,1)*69
[R14,P14,RLO14,RUP14]= corrcoef(par1, par2,'Rows','pairwise', 'alpha', 0.05)
%% 
% cor_all=horzcat([z1(2)], [z4(2)],...
%     [z7(2)],[z10(2)]);
% p_all=horzcat([p1(2) ;p2(2); p3(2)], [p4(2) ;p5(2); p6(2)],...
%     [p7(2) ;p8(2); p9(2)],[p10(2) ;p11(2); p12(2)]);
%  m=p_all<0.05;
% G=cor_all;
% % G(m==0)=m(m==0);
% % figure;imagesc(G);colorbar;
% % [cmap]=buildcmap('bwg');
% % colormap(cmap);
% % caxis([-1 1]);
%% 

%figure;bar(cor_all)
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 200, 300]);
scatter([R1(2); R4(2)],[1 3],'m','filled')
hold on;scatter([R7(2); R10(2)],[2 4],'c','filled')
ylim([0 5]);yticks([1:1:5]);xlabel('Correlation');
hold on;plot([RLO1(2) RUP1(2)],[1 1],'-m');
hold on;plot([RLO4(2) RUP4(2)],[3 3],'-m');
hold on;plot([RLO7(2) RUP7(2)],[2 2],'-c');
hold on;plot([RLO10(2) RUP10(2)],[4 4],'-c');
yticklabels({'EX','EX','IN','IN'});title('Cx L2/3 - L4');
legend('10-60°','100-150');legend boxoff;
set(gca,'FontSize',10);
%% 
