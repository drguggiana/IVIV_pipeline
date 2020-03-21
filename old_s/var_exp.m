%Display variance explained 
function var_exp(explained1,explained2,des)  
%des=name for legends
if isempty(explained2)==1 & isempty(des)==1
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y1=explained1;yyaxis left;h=bar(y1);h.EdgeColor = 'k';h.FaceColor = [1 1 1];
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:5:25]);
y3=cumsum(y1);yyaxis right
%cdfplot(r.explained);hold on;
p1=plot(y3,'--');hold on;
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[0 0 0];
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]); 
else
fig1=figure();left_color = [0 0 0];right_color = [0 0 0];set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
hold on;y2=explained1;yyaxis left;h2=bar(y2);h2.EdgeColor = 'k';h2.FaceColor = [0 0 0.7];%ylim([0 20])
hold on;y1=explained2;yyaxis left;h=bar(y1);h.EdgeColor = 'k';h.FaceColor = [0.7 0 0];
xlim([0 10.5]);ylabel('Variance explained per component');xlabel('Principal Component');
yticks([0:5:100]);
y3=cumsum(y1);y4=cumsum(y2);yyaxis right
%cdfplot(r.explained);hold on;
p1=plot(y3,'--');hold on;p2=plot(y4,'--')
xlim([0 10.5]);ylim([0 100]);ylabel('Cumulative');xlabel('Principal Component');box off;grid off;
title(''); set(gcf,'color','w');
p1.Color=[1 0 0];p2.Color=[0 0 1];legend(des{1}, des{2},des{2}, des{1});
set(gca,'FontSize',14);
yticks([0:25:100]);
xticks([1:1:10]);
end
end