function savepdf_SW(fn)
%save all open figures as PDF in a folder specified by fn
cd(fn);
 figHandles = findall(0,'Type','figure')
 for i=1:numel(figHandles)
  print(figHandles(i),'-dpdf',num2str(figHandles(i).Number));
 end
end
 
 %  for i = 1:numel(figHandles)
%      export_fig(fn, '-pdf', figHandles(i), '-append')
%  end

%  for i = 1
%      export_fig(fn, '-png', figHandles(i), '-append')
%  end
% 
%  