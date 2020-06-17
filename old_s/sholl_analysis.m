function [max_s dis_peak max_s_ba dis_peak_ba]=sholl_analysis(str,ids,pl) 

if pl==0
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

 for i=1:length(ids) 
tmp=str{1,ids(i)}  
if iscell(tmp)==1
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(tmp{1,2}, 20, '-s');
max_s_ba(i)=max(s);
temp=dd(find(s==max(s)));
dis_peak_ba(i)=temp(1);
else  
 max_s_ba(i)=NaN;   
 dis_peak_ba(i)=NaN;
end
 end
 
 
else pl==1

hold on;
for i=1:length(ids) 
tmp=str{1,ids(i)}  
if iscell(tmp)==1
fig3=figure(3)
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(tmp{1,1}, 20, '-s');
close(fig3);
max_s(i)=max(s);
temp=dd(find(s==max(s)));
dis_peak(i)=temp(1);
;p1=plot(dd,s,'-k');
p1.Color(4) = 0.1;
hold on;
else  
 max_s(i)=NaN;   
 dis_peak(i)=NaN;
end
end
hold on;

 for i=1:length(ids) 
tmp=str{1,ids(i)}  
if iscell(tmp)==1
    fig3=figure(3)
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree(tmp{1,2}, 20, '-s');
close(fig3);
max_s_ba(i)=max(s);
temp=dd(find(s==max(s)));
dis_peak_ba(i)=temp(1);
;p2=plot(dd,s,'-m');
p2.Color(4) = 0.1;
hold on;
else  
 max_s_ba(i)=NaN;   
 dis_peak_ba(i)=NaN;
end
 end
end
end