function [eye1_contra eye2_ipsi stimvideo i1 i2 v1] = getvideos(FileID,diri,what)

if nargin<2
    diri=cd;
end

if nargin<3
    what=1;
end

cd(diri)

dat = diri(end-9:end);


auxdir=['..\..\Data\' dat '\'];

if ~isnumeric(FileID)
    eye1_contra_f = dir([auxdir '\*' FileID '*.eye1*']);
    eye2_ipsi_f = dir([auxdir '\*' FileID '*.eye2*']);
    stimvideo_f = dir([auxdir '\*' FileID '*.vid*']);
else
    eye1_contra_f = dir([auxdir '\*' num2str(FileID) '*.eye1*']);
    eye2_ipsi_f = dir([auxdir '\*' num2str(FileID) '*.eye2*']);
    stimvideo_f = dir([auxdir '\*' num2str(FileID) '*.vid*']);
end

if isempty(eye1_contra_f)
    disp('LOADED PLUS ONE!')
    eye1_contra_f = dir([auxdir '\*' num2str(str2num(FileID) +1) '*.eye1*']);
    eye2_ipsi_f = dir([auxdir '\*' num2str(str2num(FileID) +1) '*.eye2*']);
    stimvideo_f = dir([auxdir '\*' num2str(str2num(FileID) +1) '*.vid*']);
end

if what > 1;
    try
        [eye1_contra i1] = load_eye_monitor_data([auxdir eye1_contra_f.name],1);
        eye1_contra = permute(eye1_contra,[2 1 3]);
    catch
        disp('No movie for eye1')
        eye1_contra=[NaN];
        i1 = [NaN];
    end
    try
        [eye2_ipsi i2] =  load_eye_monitor_data([auxdir eye2_ipsi_f.name],1);
        eye2_ipsi =  permute(eye2_ipsi,[2 1 3]);
    catch
        disp('No movie for eye2')
        eye2_ipsi=[NaN];
        i2 = [NaN];
    end
    
    if what == 3;
        try
            [stimvideo v1] =  load_vid_data([auxdir stimvideo_f.name]);
            stimvideo =  permute(stimvideo,[2 1 3]);
        catch
            disp('No stim movie / stim movie corrupted')
            stimvideo=[NaN];
            v1 = [NaN];
        end
    else
        stimvideo=[NaN];
        v1 = [NaN];
    end
    try
        %         ha = view_tiff_stack_tr(eye1_contra);
        %         set(ha, 'Name', ['ipsi: ' eye1_contra_f.name], 'NumberTitle', 'off');
    end
    try
        %         hb = view_tiff_stack_tr(eye2_ipsi);
        %         set(hb, 'Name', ['contra: ' eye2_ipsi_f.name], 'NumberTitle', 'off');
    end
    try
        %         hc = view_tiff_stack_tr(stimvideo);
        %         set(hc, 'Name', stimvideo_f.name, 'NumberTitle', 'off');
    end
    
    
elseif what ==1;
    try
        [stimvideo v1] =  load_vid_data([auxdir stimvideo_f.name]);
        eye1_contra=[NaN];
        i1 = [NaN];
        eye2_ipsi=[NaN];
        i2 = [NaN];
        %         stimvideo =  permute(stimvideo,[2 1 3]);
    catch
        disp('No stim movie')
        stimvideo=[NaN];
        v1 = [NaN];
    end
    
    try
        %         hf = view_tiff_stack_tr(stimvideo);
        %         set(hf, 'Name', stimvideo_f.name, 'NumberTitle', 'off');
        %         pause(2);
        %         truesize([size(stimvideo,1)*3 size(stimvideo,2)*3]);
    end
    
end
% view_tiff_stack_plot_aux(data, [auxdata(7,frametimes)' auxdata(8,frametimes)']);



