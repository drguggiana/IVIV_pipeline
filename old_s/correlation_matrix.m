function [G adj_p]=correlation_matrix(data,bon)

%[R,P]=corrcoef(data,'rows','complete');
%[R,P]=corrcoef(data,'rows','pairwise');
%[R,P]=corr(data,'rows','pairwise','Type','Kendall');
[R,P]=corr(data,'Type','Spearman','rows','pairwise');
%[R,P]=corrcoef(data);
 G=R;
 if bon==1
 [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P,0.05,'pdep');
 m=adj_p<0.05;
 else bon==0
  m=P<0.05;
  adj_p=m;
 end
G(m==0)=m(m==0);
G=tril(G);
figure;imagesc(G);colorbar;
[cmap]=buildcmap('bwg');
colormap(cmap);
caxis([-1 1]);
xlabel('Feature number');
ylabel('Feature number');
axis square;
set(gcf,'color','w');
set(gca,'box','off');
end