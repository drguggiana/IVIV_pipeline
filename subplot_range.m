function subplot_range(rows,columns,tar_rows,tar_cols)

%generate the vector for subplot

%allocate memory for the vector
tar_vector = zeros(length(tar_rows)*length(tar_cols),1);
%initialize a counter for the individual panels
panel_counter = 1;
%and one for the actual vector
vec_counter = 1;
%for all the rows
for row = 1:rows
    %for all the columns
    for cols = 1:columns
        %if present in the desired coordinates, include the position in the
        %vector
        if any(tar_rows==row)&&any(tar_cols==cols)
            %save the current panel position
            tar_vector(vec_counter) = panel_counter;
            %update the vector counter
            vec_counter =vec_counter + 1;
        end
        %update the panel counter
        panel_counter = panel_counter + 1;
    end
end

%call subplot with the corresponding arguments
subplot(rows,columns,tar_vector)