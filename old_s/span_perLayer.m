function [ex_spanhL23,ex_spanhL4,ex_spanhL5,in_spanhL23,in_spanhL4,in_spanhL5] = span_perLayer(ex_map,in_map,nan_vector)
incl_idx=1;
for i=1:length(nan_vector(incl_idx:end))
    tmp1=ex_map(3:5,:,i);
    for k=1:16
        tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
    end
    tmp2=in_map(3:5,:,i);
    for k=1:16
        tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
    end
    tmp=find(tL23e>0);
    ex_spanhL23(i)=tmp(end)-tmp(1);
    tmp=[];
    tmp=find(tL23i>0);
    in_spanhL23(i)=tmp(end)-tmp(1);
    tmp=[];
    tL23i=[];
    tL23e=[];
    tmp1=[];tmp2=[];
end
for i=1:length(nan_vector(incl_idx:end))
    tmp1=ex_map(6:7,:,i);
    for k=1:16
        tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
    end
    tmp2=in_map(6:7,:,i);
    for k=1:16
        tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
    end
    tmp=find(tL23e>0);
    if isempty(tmp)==0
        ex_spanhL4(i)=tmp(end)-tmp(1);
    else
        ex_spanhL4(i)=NaN;
    end
    tmp=[];
    tmp=find(tL23i>0);
    if isempty(tmp)==0
        in_spanhL4(i)=tmp(end)-tmp(1);
    else
        in_spanhL4(i)=NaN;
    end
    tmp=[];
    tL23i=[];
    tL23e=[];
end
%L5
for i=1:length(nan_vector(incl_idx:end))
    tmp1=ex_map(8:10,:,i);
    for k=1:16
        tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
    end
    tmp2=in_map(8:10,:,i);
    for k=1:16
        tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
    end
    tmp=find(tL23e>0);
    if isempty(tmp)==0
        ex_spanhL5(i)=tmp(end)-tmp(1);
    else
        ex_spanhL5(i)=NaN;
    end
    tmp=[];
    tmp=find(tL23i>0);
    if isempty(tmp)==0
        in_spanhL5(i)=tmp(end)-tmp(1);
    else
        in_spanhL5(i)=NaN;
    end
    tmp=[];
    tL23i=[];
    tL23e=[];
end
end