function [frac_exh,abs_exh,frac_inh,abs_inh,frac_exv,abs_exv,frac_inv,abs_inv,L23h,L4h,L5h,pialD,layers_assign] = iviv_profiles(u,str) 

%Inputs 
%u= idx of ipsi, contra, bino or unresponsive cells 
%str= saved structure on R that contains iviv information 

%Outputs
%fraction of excitation/inhbition per row/column
%absolute sum of excitation/inhbition per row/column

% allocate memory
exmap = zeros(16,16,length(u));
inhmap = zeros(16,16,length(u));
layers = zeros(16,16,length(u));

frac_exh = zeros(length(u),16);
abs_exh = zeros(length(u),16);

frac_inh = zeros(length(u),16);
abs_inh = zeros(length(u),16);

frac_exv = zeros(length(u),16);
abs_exv = zeros(length(u),16);

frac_inv = zeros(length(u),16);
abs_inv = zeros(length(u),16);

layers_assign = zeros(14,length(u));
%Read out 16 x 16 input maps per cell 
for i=1:length(u)
    
    if ~isempty(str(u(i)).subpixel_excMap)
        exmap(:,:,i)=str(u(i)).subpixel_excMap;
    else
        exmap(:,:,i)=NaN;
    end
    if ~isempty(str(u(i)).subpixel_excMap)
        inhmap(:,:,i)=str(u(i)).subpixel_inhMap;
    else
        inhmap(:,:,i)=NaN;
    end
    if size(str(u(i)).layers) == [16 16]
        layers(:,:,i)=str(u(i)).layers;
    else
        layers(:,:,i)=NaN;
    end
end

%HORIZONTAL Fraction excitation
for i=1:length(u)
    if isempty(str(u(i)).subpixel_excMap)
        frac_exh(i,1:16)=NaN;
        abs_exh(i,1:16)=NaN;
        continue
    end
    %only use 2:14 rows cause of L1 and WM exclusion
    for k=1:16%for horizontal input include all
        temp=exmap(3:16,:,i) ;
        %layers_assign(:,1,i)=layers(3:16,1,i);
        frac_exh(i,k)=sum(temp(:,k))/sum(temp(:));
        abs_exh(i,k)=sum(temp(:,k));
    end
end

%HORIZONTAL Fraction Inhibition
for i=1:length(u)
    if isempty(str(u(i)).subpixel_inhMap)
        frac_inh(i,1:16)=NaN;
        abs_inh(i,1:16)=NaN;
        continue
    end
    for k=1:16
        temp=inhmap(1:16,:,i) ;
        frac_inh(i,k)=sum(temp(:,k))/sum(temp(:));
        abs_inh(i,k)=sum(temp(:,k));
    end
end


%VERTICAL Fraction excitation
for i=1:length(u)
    if isempty(str(u(i)).subpixel_excMap)
        layers_assign(:,i)=NaN;
        frac_exv(i,1:16)=NaN;
        abs_exv(i,1:16)=NaN;
        continue
    end
    for k=1:14
        temp=exmap(3:16,:,i) ;
        layers_assign(:,i)=layers(3:16,1,i);
        frac_exv(i,k)=sum(temp(k,:))/sum(temp(:));
        abs_exv(i,k)=sum(temp(k,:));
    end
end

%VERTICAL Fraction inhibition
for i=1:length(u)
    if isempty(str(u(i)).subpixel_inhMap)
        frac_inv(i,1:16)=NaN;
        abs_inv(i,1:16)=NaN;
        continue
    end
    for k=1:16
        temp=inhmap(1:16,:,i) ;
        %layers_assign(:,1,i)=layers(1:16,1,i);
        frac_inv(i,k)=sum(temp(k,:))/sum(temp(:));
        abs_inv(i,k)=sum(temp(k,:));
    end
end

pialD=[str(u).pialD]';

%% Horizontal fraction for L23 L4 and L5

% for i=1:length(u)
%     if ~isempty(str(u(i)).subpixel_excMap)
%         %L23ex
%         temp=exmap(3:5,:,i);
%         L23exh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%         %L4ex
%         temp=exmap(6:7,:,i);
%         L4exh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%         %L5ex
%         temp=exmap(8:10,:,i);
%         L5exh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%     end
%     if ~isempty(str(u(i)).subpixel_inhMap)
%         %L23in
%         temp=inhmap(3:5,:,i);
%         L23inh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%         %L4in
%         temp=inhmap(6:7,:,i);
%         L4inh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%         %L5in
%         temp=inhmap(8:10,:,i);
%         L5inh(i,:)=[sum(sum(temp(:,1:8)))/sum(temp(:)) sum(sum(temp(:,9:16)))/sum(temp(:))];
%     end
% end
% L23h=[L23exh L23inh];
% L4h=[L4exh L4inh];
% L5h=[L5exh L5inh];

L23h=[];
L4h=[];
L5h=[];
end