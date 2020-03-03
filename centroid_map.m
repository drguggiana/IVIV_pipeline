function [out_ang] = centroid_map(test_map,somax,pia_input,idx,row_shift);
%Calculate weighted centroid and total mass
for i=1:length(test_map)
A=test_map(:,:,i);
tot_mass = sum(A(:));
[ii,jj] = ndgrid(1:size(A,1),1:size(A,2));
  R = sum(ii(:).*A(:))/tot_mass;
  C = sum(jj(:).*A(:))/tot_mass;
%props2 = regionprops(true(size(test_map(:,:,i))), test_map(:,:,i), 'WeightedCentroid');
% props2=regionprops(A,'Centroid');
% if isempty(props2)==0;
%  C=props2.Centroid(1);
%  R=props2.Centroid(2);
% else
%    C=NaN;
%    R=NaN;
% end
out(i,:) = [tot_mass,R,C];
cell_cord(i,:)=[(8+somax(idx(i))/69) pia_input(idx(i))/69];
%val_hot(i)=A(round(R),round(C));
end
%ALTERNATIVE for centroid
%props2 = regionprops(true(size(test_map(:,:,i))), test_map(:,:,i), 'WeightedCentroid');
% x_weight(i)=props2.WeightedCentroid(1);
% y_weight(i)=props2.WeightedCentroid(2);

wx=out(:,3);
wy=out(:,2)+row_shift;

sx=cell_cord(:,1);
sy=cell_cord(:,2);
dpx=wx-sx;
dpy=wy-sy;
%Calculate alpha and beta
for i=1:length(test_map)
hypo(i)=sqrt(((dpx(i)^2)+(dpy(i)^2)))
sina(i)=dpy(i)/hypo(i);
ang_a(i)=asind(sina(i));
ang_b(i)=90-ang_a(i);
centr_vector=polyfit([sx(i), wx(i)],[sy(i), wy(i)],1);
vec_slope(i)=centr_vector(1);
if dpy(i)<0 & dpx(i)<0 %top left
ang_b_v(i)=ang_b(i)*-1;
ang_a_v(i)=ang_a(i)*-1;
qd(i)=1;
elseif dpy(i)>0 & dpx(i)<0 %bottom left
ang_b_v(i)=ang_b(i)*-1;  
ang_a_v(i)=ang_a(i)*-1;
qd(i)=2;
elseif dpy(i)>0 & dpx(i)>0 %bottom right
ang_b_v(i)=ang_b(i);
ang_a_v(i)=ang_a(i);
qd(i)=3;
else dpy(i)<0 & dpx(i)>0 %top right
ang_b_v(i)=ang_b(i); 
ang_a_v(i)=ang_a(i);
qd(i)=4;
end
end
%out
out_ang=[sx sy wx wy ang_a' ang_b' ang_a_v' hypo' vec_slope' qd'];
end