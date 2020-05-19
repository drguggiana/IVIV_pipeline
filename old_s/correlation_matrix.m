function G=correlation_matrix(data,bon)

%[R,P]=corrcoef(data);
[R,P]=corrcoef(data,'Rows','pairwise')
 G=R;
 if bon==1
 [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P,0.05,'pdep');
 m=adj_p<crit_p;
 else bon==0
  m=P<0.05;
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