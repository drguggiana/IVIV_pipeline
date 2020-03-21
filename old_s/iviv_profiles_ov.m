function [frac_ovh abs_ovh frac_ovv abs_ovv] = iviv_profiles_ov(u);
%Inputs 


%Outputs
%fraction of excitation/inhbition per row/column
%absolute sum of excitation/inhbition per row/column
%HORIZONTAL Fraction 

for i=1:length(u);
%only use 2:14 rows cause of L1 and WM exclusion 

for k=1:16%for horizontal input include all
temp=u(1:16,:,i) ;
%layers_assign(:,1,i)=layers(3:16,1,i);
frac_ovh(i,k)=sum(temp(:,k))/sum(abs(temp(:)));
abs_ovh(i,k)=sum(temp(:,k));
end
end
temp=[];

%VERTICAL Fraction 
frac_inv=[];
abs_inv=[];
for i=1:length(u);
for k=1:16
temp=u(1:16,:,i) ;
%layers_assign(:,1,i)=layers(1:16,1,i);
frac_ovv(i,k)=sum(temp(k,:))/sum(abs(temp(:)));
abs_ovv(i,k)=sum(temp(k,:));
end
end

%layers_ass=layers_assign(:,:)';



end 