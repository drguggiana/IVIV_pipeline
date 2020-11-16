function [L1h L23h L4h L5h] = h_fraclayer(ex_map, in_map)

map=ex_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(1:2,1:16,i);
ex(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];

map=in_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(1:2,1:16,i);
in(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];
L1h=[ex in];

map=ex_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(3:5,1:16,i);
ex(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];

map=in_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(3:5,1:16,i);
in(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];
L23h=[ex in];

map=ex_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(6:7,1:16,i);
ex(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];

map=in_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(6:7,1:16,i);
in(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];
L4h=[ex in];


map=ex_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(8:11,1:16,i);
ex(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];

map=in_map;
for i=1:length(map)
my=[];mx=[];
temp=[];
temp=map(8:11,1:16,i);
in(i,:)=sum(temp)/sum(sum(temp));
end
temp=[];
L5h=[ex in];
end