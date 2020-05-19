
function [max_p]=max_pixel(ex_map, in_map)
for i=1:length(ex_map)
[iny inx]=find(in_map(:,:,i)==max(max(in_map(:,:,i))));
[exy exx]=find(ex_map(:,:,i)==max(max(ex_map(:,:,i))));
cexx(i)=exx(1);
cexy(i)=exy(1);
cinx(i)=inx(1);
ciny(i)=iny(1);
end
max_p=[cexx ;cexy; cinx; ciny;];
end