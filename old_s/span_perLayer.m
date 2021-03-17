function [ex_spanhL23,ex_spanhL4,ex_spanhL5,in_spanhL23,in_spanhL4,in_spanhL5] = span_perLayer(str)

% get the number of cells
cell_num = length(str);


% allocate memory for the outputs
ex_spanhL23 = zeros(cell_num,1);
in_spanhL23 = zeros(cell_num,1);

% for all the cells
for i=1:cell_num
    
    tL23e = zeros(16,1);
    tL23i = zeros(16,1);
    % load the excitation map
    if ~isempty(str(i).subpixel_excMap) && ~all(isnan(str(i).subpixel_excMap(:)))
        tmp1=str(i).subpixel_excMap(3:5,:);
        for k=1:16
            tL23e(k)=sum(tmp1(:,k))/sum(tmp1(:));
        end
        tmp=find(tL23e>0);
        ex_spanhL23(i)=tmp(end)-tmp(1);
    else
        ex_spanhL23(i) = NaN;
    end
    if ~isempty(str(i).subpixel_inhMap) && ~all(isnan(str(i).subpixel_inhMap(:)))
        tmp2=str(i).subpixel_inhMap(3:5,:);
        for k=1:16
            tL23i(k)=sum(tmp2(:,k))/sum(tmp2(:));
        end
        tmp=find(tL23i>0);
        in_spanhL23(i)=tmp(end)-tmp(1);
    else
        in_spanhL23(i) = NaN;
    end
end

%L4
ex_spanhL4 = zeros(cell_num,1);
in_spanhL4 = zeros(cell_num,1);

for i=1:cell_num
    
    tL4e = zeros(16,1);
    tL4i = zeros(16,1);
    if ~isempty(str(i).subpixel_excMap) && ~all(isnan(str(i).subpixel_excMap(:)))
        tmp1=str(i).subpixel_excMap(6:7,:);
        for k=1:16
            tL4e(k)=sum(tmp1(:,k))/sum(tmp1(:));
        end
        tmp=find(tL4e>0);
        if isempty(tmp)==0
            ex_spanhL4(i)=tmp(end)-tmp(1);
        else
            ex_spanhL4(i)=NaN;
        end
    else
        ex_spanhL4(i) = NaN;
    end
    if ~isempty(str(i).subpixel_inhMap) && ~all(isnan(str(i).subpixel_inhMap(:)))
        tmp2=str(i).subpixel_inhMap(6:7,:);
        for k=1:16
            tL4i(k)=sum(tmp2(:,k))/sum(tmp2(:));
        end
        tmp=find(tL4i>0);
        if isempty(tmp)==0
            in_spanhL4(i)=tmp(end)-tmp(1);
        else
            in_spanhL4(i)=NaN;
        end
    else
        ex_spanhL4(i) = NaN;
    end
end

%L5
ex_spanhL5 = zeros(cell_num,1);
in_spanhL5 = zeros(cell_num,1);

for i=1:cell_num
    tL5e = zeros(16,1);
    tL5i = zeros(16,1);
    if ~isempty(str(i).subpixel_excMap) && ~all(isnan(str(i).subpixel_excMap(:)))
        tmp1=str(i).subpixel_excMap(8:10,:);
        for k=1:16
            tL5e(k)=sum(tmp1(:,k))/sum(tmp1(:));
        end
        tmp=find(tL5e>0);
        if isempty(tmp)==0
            ex_spanhL5(i)=tmp(end)-tmp(1);
        else
            ex_spanhL5(i)=NaN;
        end
    else
        ex_spanhL5(i) = NaN;
    end
    if ~isempty(str(i).subpixel_inhMap) && ~all(isnan(str(i).subpixel_inhMap(:)))
        tmp2=str(i).subpixel_inhMap(8:10,:);
        for k=1:16
            tL5i(k)=sum(tmp2(:,k))/sum(tmp2(:));
        end
        
        tmp=find(tL5i>0);
        if isempty(tmp)==0
            in_spanhL5(i)=tmp(end)-tmp(1);
        else
            in_spanhL5(i)=NaN;
        end
    else
        in_spanhL5(i) = NaN;
    end
end
end