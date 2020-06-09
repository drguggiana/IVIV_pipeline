function [out_ang] = centroid_map(test_map,somax,somay,pialD,row_shift,varargin)
% calculate the centroid and centroid related metrics for the given set of
% maps

% parse the varargin
if length(varargin) >= 1
    % take the given constant
    we = varargin{1};
else
    % default is weighted
    we = 1;
end

% correct the row_shift
row_shift = row_shift + 0.5;

% allocate memory
out = zeros(length(test_map),3);
cell_coord = zeros(length(test_map),2);

%Calculate weighted centroid and total mass
for i = 1:length(test_map)
    % load the current map
    A = vertcat(mean(test_map(1:2,:,i),1),test_map(2:end-1,:,i),mean(test_map(end-1:end,:,i),1));

    % calculate the total mass
    tot_mass = sum(A(:));
    % calculate the weighted centroid
    [ii,jj] = ndgrid([0,0.5:size(test_map,1)-2.5,size(test_map,1)-2],1:size(A,2));
%     [ii,jj] = meshgrid(1:size(A,1)+1,1:size(A,2)+1);
%     A = interp2(ii0,jj0,A',ii,jj);
    
    R = nansum(ii(:).*A(:))/tot_mass;
    C = nansum(jj(:).*A(:))/tot_mass;
    % store the masses in total and every dimension
    if we == 1
        out(i,:) = [tot_mass,R,C];
        
    else 
        tt = regionprops(abs(A)>0,'Centroid','Area');
        if isempty(tt)
            out(i,:) = [tot_mass, NaN,NaN];

        elseif length(tt) > 1
            [~,max_idx] = max(vertcat(tt.Area));
            out(i,:) = [tot_mass, tt(max_idx).Centroid(2),tt(max_idx).Centroid(1)];
        else
            out(i,:) = [tot_mass, tt.Centroid(2),tt.Centroid(1)];
        end

    end
    % store the soma coordinates (converting to grid pixels and 
    % reversing the y axis)
%     cell_coord(i,:)=[(8.5+somax(idx(i))/69) pia_input(idx(i))/69];
    cell_coord(i,1) = 8.5 + somax(i)/69;
    cell_coord(i,2) = 8.5 - somay(i)/69;
end

% get the centroid coordinates and convert to grid pixels (from image
% pixels)
% wx = out(:,3) - 8.5;
% % wy = out(:,2)+row_shift;
% wy = out(:,2) - 8.5 + row_shift;
wx = out(:,3);
wy = out(:,2) + row_shift;
% get the soma coordinates
sx = cell_coord(:,1);
sy = cell_coord(:,2);
% correct the y centroid and soma position
soma_delta = pialD./69 - sy;
sy = pialD./69;
wy = wy + soma_delta;

% get the delta
dpx = wx-sx;
dpy = wy-sy;

% allocate memory
hypo = zeros(length(test_map),1);
sina = zeros(length(test_map),1);
ang_a = zeros(length(test_map),1);
ang_b = zeros(length(test_map),1);
vec_slope = zeros(length(test_map),1);

ang_b_v = zeros(length(test_map),1);
ang_a_v = zeros(length(test_map),1);
qd = zeros(length(test_map),1);


%Calculate alpha and beta
for i=1:length(test_map)
    hypo(i)=sqrt(((dpx(i)^2)+(dpy(i)^2)));
    sina(i)=dpy(i)/hypo(i);
%     ang_a(i)=asind(sina(i));
    ang_a(i) = rad2deg(atan2(dpy(i),dpx(i)));
    ang_b(i)=90-ang_a(i);
    centr_vector=polyfit([sx(i), wx(i)],[sy(i), wy(i)],1);
    vec_slope(i)=centr_vector(1);
    if dpy(i)<0 && dpx(i)<0 %top left
        ang_b_v(i)=ang_b(i)*-1;
        ang_a_v(i)=ang_a(i)*-1;
        qd(i)=1;
    elseif dpy(i)>0 && dpx(i)<0 %bottom left
        ang_b_v(i)=ang_b(i)*-1;
        ang_a_v(i)=ang_a(i)*-1;
        qd(i)=3;
    elseif dpy(i)>0 && dpx(i)>0 %bottom right
        ang_b_v(i)=ang_b(i);
        ang_a_v(i)=ang_a(i);
        qd(i)=4;
    else %dpy(i)<0 && dpx(i)>0 %top right
        ang_b_v(i)=ang_b(i);
        ang_a_v(i)=ang_a(i);
        qd(i)=2;
    end
end
%output
out_ang=[sx sy wx wy ang_a ang_b ang_a_v hypo vec_slope qd];