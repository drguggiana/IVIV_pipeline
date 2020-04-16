function [frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv L23h L4h L5h pialD layers_assign] = iviv_profiles(u,ex_map,in_map) 

%Read out 16 x 16 input maps per cell 
for i=1:length(u);
exmap(:,:,i)=ex_map(:,:,i);
inhmap(:,:,i)=in_map(:,:,i);

end
%HORIZONTAL Fraction excitation
frac_exh=[];
abs_exh=[];
for i=1:length(u);
%only use 2:14 rows cause of L1 and WM exclusion 
for k=1:16%for horizontal input include all
temp=exmap(3:16,:,i) ;
%layers_assign(:,1,i)=layers(3:16,1,i);
frac_exh(i,k)=sum(temp(:,k))/sum(temp(:));
abs_exh(i,k)=sum(temp(:,k));
end
end
temp=[];
%HORIZONTAL Fraction Inhibition 
frac_inh=[];
abs_inh=[];
for i=1:length(u);
for k=1:16
temp=inhmap(1:16,:,i) ; 
frac_inh(i,k)=sum(temp(:,k))/sum(temp(:));
abs_inh(i,k)=sum(temp(:,k));
end
end


%VERTICAL Fraction excitation
frac_exv=[];
abs_exv=[];
for i=1:length(u);
for k=1:14
temp=exmap(3:16,:,i) ;
%layers_assign(:,i)=layers(3:16,1,i);
frac_exv(i,k)=sum(temp(k,:))/sum(temp(:));
abs_exv(i,k)=sum(temp(k,:));
end
end
%VERTICAL Fraction inhibition
frac_inv=[];
abs_inv=[];
for i=1:length(u);
for k=1:16
temp=inhmap(1:16,:,i) ;
%layers_assign(:,1,i)=layers(1:16,1,i);
frac_inv(i,k)=sum(temp(k,:))/sum(temp(:));
abs_inv(i,k)=sum(temp(k,:));
end
end
end