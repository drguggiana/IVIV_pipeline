function [idata, imeta_info] = load_eye_monitor_data(file_name, load_img)


if nargin < 1
    disp('usage: load_labview_acute(file_name[, file_format])');
    return;
end

meta_info_size=9;

%% try to open file, otherwise quit
[fi, message] = fopen(file_name, 'r', 'ieee-be');

status = 0;
data = [];
scanrateA = [];

if fi == -1
    disp('There was a problem reading the following file:');
    disp(file_name);
    disp(message);
    status = -1;
    return;
end

%% get file size

meta_info = fread(fi, meta_info_size, 'double');
if meta_info(3)==0 | rem(meta_info(3),1)~=0
    size_x=meta_info(5);
    size_y=meta_info(6);
else
    size_x=meta_info(5)-meta_info(3);
    size_y=meta_info(6)-meta_info(4);
end
fseek(fi, 0, 'eof');
file_size = ftell(fi);
fclose(fi);


%% open the file again - this time to read in the data
fi = fopen(file_name, 'r', 'ieee-be');

%% get all header information


% double = 8 bytes, uint8 = 1 byte
nbr_frames=(file_size)/(size_x*size_y+meta_info_size*8);
if rem(nbr_frames,1) == 0 ; % check if integer frame number - otherwise the movie is corrupted!
    if load_img
        idata=zeros(size_x,size_y,nbr_frames,'uint8');
    else
        idata=[];
    end
else disp('Non-integer number of frames detected. Movie is most likely corrupted!')
    fclose(fi);
    return
end

imeta_info=zeros(9,nbr_frames);

fprintf(1,'loading eyemovie frame (of %d)     ', nbr_frames);
go_on=1;
for ind=1:nbr_frames
    fprintf(1,'\b\b\b\b\b%5.5s',num2str(ind));
    imeta_info(:,ind) = fread(fi, meta_info_size, 'double');
    if load_img
        idata(:,:,ind) = uint8(reshape(fread(fi, size_x*size_y, 'uint8'),size_x,size_y));
    else
        fseek(fi,size_x*size_y,'cof');
    end
end
fprintf('\n')


%% read in the data

fclose(fi);


%% EOF
