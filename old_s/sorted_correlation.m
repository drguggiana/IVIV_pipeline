function sorted_correlation(data,sort_data)

com=[];com=[data sort_data];
[R1,P1,RLO1,RUP1]=corrcoef(com,'rows','pairwise');
corr_tree=[];k=[];
[corr_tree k]=sort(R1(end,1:end-1));
temp1=[];temp2=[];temp3=[];
temp1=RLO1(end,1:end-1);
temp2=RUP1(end,1:end-1);
up_tree=[];low_tree=[];
up_tree=temp1(k);
low_tree=temp2(k);
temp3=[];temp3=P1(end,1:end-1);
p_tree=temp3(k);

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 200, 400])
for i=1:length(corr_tree)
hold on;scatter(corr_tree(i),i,'m','filled')
hold on;plot([low_tree(i) up_tree(i)],[i i],'-m');
if p_tree(i)<0.001
    text(1.05,i,'***')
elseif p_tree(i)<0.01
    text(1.05,i,'**')
elseif p_tree(i)<0.05
     text(1.05,i,'*')
else
end
end
%ylim([0 10]);

hold on;line([0 0], [0 10],'Color','k','LineStyle','--');
% xlabel('Correlation with pial depth');title('Apical tree');
% stri={'RDA_{max}(µm)','Total Length (µm)','Peak Nr. crossing','Dis. peak branch (µm)','WHA',...
%      'RDB_{max}(µm)','Total Length (µm)','Peak Nr. crossing','Dis. peak branch(µm)','WHA'}
end