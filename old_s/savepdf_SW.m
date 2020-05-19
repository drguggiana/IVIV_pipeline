function savepdf_SW(fn,form)
%save all open figures as PDF in a folder specified by fn
cd(fn);
 figHandles = findall(0,'Type','figure')
 if form==1
 for i=1:numel(figHandles)
  print(figHandles(i),'-dpdf',num2str(figHandles(i).Number));
 
 end
 else
    for i=1:numel(figHandles)
 savefig(figHandles(i),num2str(figHandles(i).Number));
 
 end
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