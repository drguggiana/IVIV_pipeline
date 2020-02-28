function correlation_matrix(data,bon)

[R,P]=corrcoef(data,'rows','complete');
%[R,P]=corrcoef(data);
 G=R;
 if bon==1
 [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P,0.05,'pdep');
 m=adj_p<0.05;
 else bon==0
  m=P<0.05;
 end
G(m==0)=m(m==0);
G=tril(G);
figure;imagesc(G);colorbar;
[cmap]=buildcmap('bwr');
colormap(cmap);
caxis([-1 1]);
xlabel('Feature number');
ylabel('Feature number');
axis square;
set(gcf,'color','w');
set(gca,'box','off');
end