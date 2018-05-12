function mat_out = normr_2(mat_in,varargin)

%function to normalize matrices, depending on the selection of flag.
%Omitting the flag normalizes the entire matrix

%if there are no extra arguments
if nargin == 1
    %default to normalizing the whole matrix
    flag_num = 0;
else %otherwise
    %assign whatever flag the user selected
    flag_num = varargin{1};
end
%switch normalizing depending on the user selection
switch flag_num
    %case 0 and default, normalize the whole matrix
    case 0
        
        mat_out = (mat_in - min(mat_in(:)))./(max(mat_in(:)) - min(mat_in(:)));
    %normalize each row
    case 1
        %allocate memory for the output
        mat_out = zeros(size(mat_in));
        %for all the rows
        for rows = 1:size(mat_in,1)
            
            mat_out(rows,:) = (mat_in(rows,:) - min(mat_in(rows,:)))/(max(mat_in(rows,:)) - min(mat_in(rows,:)));
        end
    %normalize each column
    case 2
        %allocate memory for the output
        mat_out = zeros(size(mat_in));
        %for all the rows
        for cols = 1:size(mat_in,2)
            
            mat_out(:,cols) = (mat_in(:,cols) - min(mat_in(:,cols)))/(max(mat_in(:,cols)) - min(mat_in(:,cols)));
        end
end
