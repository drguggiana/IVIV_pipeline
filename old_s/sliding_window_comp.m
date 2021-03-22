function [p1 p2]=sliding_window_comp(in_vec,comp_vec)

%% normalize the input vector
tr = rescale(in_vec);
[y k] = sort(tr);
%take the lower and upper 25%
usec=find(y<0.25);
lsec=find(y>0.75);
%% sort the comp_vec accordingly to the in_vec       
y2=comp_vec(k);
%sliding window size (10%)
samplesCounter = round(length(tr)*0.1);
indexes = 1 : samplesCounter : length(y);
for k = 2 : length(indexes)
	index1 = indexes(k-1);
	index2 = indexes(k) - 1;
	% Option 1 : Everything in between, inclusive.
	theSum1(k) = mean(y2(index1:index2)); 
    valu(:,k)=y2(index1:index2);
	% Option 2 : Sum y only at those two indexes, and not everything in between.
% 	theSum2(k) = 0.5 * y2(index1) + 0.5 * y2(index2);
end

%% Compare groups against upper and lower quartal
for i=1:size(valu,2)
p1(i) = ranksum(y2(usec),valu(:,i))
p2(i) = ranksum(y2(lsec),valu(:,i))
end
figure;set(gcf,'color','w');plot(1:size(valu,2),p1,'r--*');
hold on;box off;
plot(1:size(valu,2),p2,'b--*');xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1' })
hold on;line([1 size(valu,2)], [0.05 0.05],'Color','k','LineStyle','--');
ylabel('p values');
xlabel('relative pial depth');
legend('upper comp','lower comp');
legend boxoff
end