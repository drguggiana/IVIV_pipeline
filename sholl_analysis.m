function [max_s dis_peak]=sholl_analysis(str,ids) 
for i=1:length(ids) 
tmp=str{1,ids(i)}  
if iscell(tmp)==1
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(tmp{1,1}, 20, '-s');
max_s(i)=max(s);
temp=dd(find(s==max(s)));
dis_peak(i)=temp(1);
else  
 max_s(i)=NaN;   
 dis_peak(i)=NaN;
end
end
end