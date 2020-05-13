function mask = create_mask(X,Y,center_x,center_y,radius)
% create a circular mask of the given radius and center

mask = (Y - center_y).^2 + (X - center_x).^2 <= radius.^2;