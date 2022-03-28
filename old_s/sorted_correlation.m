function [corr_tree p_tree]=sorted_correlation(data,sort_data,stri,col,bon)

com=[];com=[data sort_data];
%[R1,P1,RLO1,RUP1]=corrcoef(com,'rows','pairwise');
[R1,P1,RLO1,RUP1]=corr(com,'Type','Spearman','Rows','pairwise');
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
lab_tree=stri(k);

  if bon==1
  [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_tree,.05,'pdep','yes');
% %  [adj_p, c_alpha, h] = fdr_BH(p_tree, 0.05, false); 
% %  %[adj_p_bon, c_alpha, h]=fwer_bonf(p_tree, 0.05, false);
% %  %  [adj_p, c_alpha, h]= fwer_sidak(p_tree, 0.05, false);
  p_tree=adj_p;
  else bon==0
  p_tree=p_tree;
  end

fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 400, 600])
for i=1:length(corr_tree)
hold on;scatter(corr_tree(i),i,'MarkerFaceColor',col,'MarkerEdgeColor',col)
hold on;plot([low_tree(i) up_tree(i)],[i i],'Color',col);

%  if bon==1
%   [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_tree,0.05,'pdep');
%   p_tree=adj_p;
%  else bon==0
%  p_tree=p_tree;
%   end
%  p_all(:,i)= p_tree;
 
if p_tree(i)<0.001
    text(1.05,i,'***')
elseif p_tree(i)<0.01
    text(1.05,i,'**')
elseif p_tree(i)<0.05
     text(1.05,i,'*')
else
end

end

xlim([-1 1.1]);
ylim([0 length(p_tree)]);
yticklabels([' ' lab_tree]);
yticks([0:1:length(p_tree)]);
hold on;line([0 0], [0 length(corr_tree)],'Color','k','LineStyle','--');

% xlabel('Correlation with pial depth');title('Apical tree');
% stri={'RDA_{max}(�m)','Total Length (�m)','Peak Nr. crossing','Dis. peak branch (�m)','WHA',...
%      'RDB_{max}(�m)','Total Length (�m)','Peak Nr. crossing','Dis. peak branch(�m)','WHA'}
end