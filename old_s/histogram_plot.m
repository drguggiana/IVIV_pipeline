
% 
% %histogram plot for all parameters independent of cluster yet
% com=[w.com;w2.com;w3.com];
% pia=[kk(:,1);mm(:,1);nn(:,1)];
% com=[com pia];
% 
% traces=[w.traces;w2.traces;w3.traces];
% traces=[traces num2cell(pia) num2cell(ID)];
% traces=[traces num2cell(r.idx)];

figure();
 str={'RDA_{max}(µm)','LA_{total} (µm)','PLA_{max} (µm)','BPA','BOA_{max}','BLA (µm)','PLA (µm)','WHA','WZA','XSA (µm)','YSA (µm)','ZSA (µm)',...
     'RDB_{max}(µm)','LB_{total} (µm)','PLB_{max} (µm)','BPB','BOB_{max}','BLB (µm)','PLB (µm)','WHB','WZB','XSB (µm)','YSB (µm)','ZSB (µm)', 'NB','Pial depth (µm)'}
for i=1:26
hold on;
subplot(5,6,i)
h=histogram(com(:,i),8,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5)
h.EdgeColor = 'k';
h.FaceColor = [0.5 0.5 0.5];
xlabel(str(i));
ylim([0 150]);
xlim([0 max(com(:,i))+max(com(:,i))*0.25]);
hAxis = gca;
hAxis.YAxisLocation = 'left';    % 'left' (default) or 'right'
hAxis.XAxisLocation = 'bottom'
box off
end



[R,P]=corrcoef(com,'rows','complete');
% m=P<0.05;
 G=R;
 [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P);
 m=adj_p<crit_p;
G(m==0)=m(m==0);
G=tril(G);
figure;imagesc(G);colorbar;
colormap(cool);
caxis([-1 1]);
xlabel('Feature number');
ylabel('Feature number');
set(gca, 'XTick', [1:26]);
set(gca, 'YTick', [1:26]);

textStrings = num2str(G(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:length(G));  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(G(:) > midValue, 1, 3);  