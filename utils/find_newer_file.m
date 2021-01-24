function newest_path = find_newer_file(target_folder)
% find the newest file in the target folder

% get a list of the files
file_list = dir(target_folder);
file_list = file_list([file_list.isdir]==0);

% get the number of files
file_num = length(file_list);
% allocate memory for the file dates
file_times = zeros(file_num,1);

% for all the files
for files = 1:file_num
    
    % get the file dates
    temp = GetFileTime(fullfile(file_list(files).folder,file_list(files).name));
    % get only the creation date
    file_times(files) = datenum(temp.Creation);
end

% select the newest one and output that path
[~,target_idx] = max(file_times);

newest_path = fullfile(file_list(target_idx).folder,file_list(files).name);