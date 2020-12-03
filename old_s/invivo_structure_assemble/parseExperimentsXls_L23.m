function [batchopt] = parseExperimentsXls_L23(path)

[xls_num,xls_txt]=xlsread(path);

loadcol        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'BatchAnalyze')));
mousecol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Animal_ID')));
experiment     = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExperimentalDay')));
loadrivecol      =find(~cellfun(@isempty, strfind(xls_txt(1,:),'loaddrive')));
binocol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExpIDbinoc')));
sponcol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExpIDsSpont')));
sftfcol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExpIDsSFTF')));
invitro_rec   = find(~cellfun(@isempty, strfind(xls_txt(1,:),'invitroRecordings')));
ivivROI       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ivivROI')));
invitrob    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'invitro')));
cellid      = find(~cellfun(@isempty, strfind(xls_txt(1,:),'cellID')));




k = 1;

batchopt.XLS.txt = xls_txt;
batchopt.XLS.num = xls_txt;

for i = 2:size(xls_txt,1)
    ana{k}= xls_num(i-1,loadcol-2);
    
    
    
    if ~ana{k}
        disp(['skipping experiments ' xls_txt{i,mousecol} ' (no batchload flag)']);
        continue
    end
    
    batchopt.mouse{k} = xls_txt(i,mousecol);
    batchopt.exp_invitro{k} = xls_txt(i,experiment);
   
    
    expcellids{k}                = xls_txt(i,binocol );
    expcellids2{k}               = xls_txt(i,sponcol);
    expcellids3{k}               = xls_txt(i,sftfcol);
    expcellids4{k}               = xls_txt(i,invitro_rec);
    expcellids5{k}               = xls_txt(i,ivivROI);
    expcellids6{k}               = xls_txt(i,cellid);
     
    
    batchopt.binoexp_ids{k}         = str2num((expcellids{k}{1}));
    batchopt.sponexp_ids{k}         = str2num((expcellids2{k}{1}));
    batchopt.sftfexp_ids{k}         = str2num((expcellids3{k}{1}));
    batchopt.invitro_rec{k}          = str2num((expcellids4{k}{1}));
    batchopt.ivivROI{k}             = str2num((expcellids5{k}{1}));
    batchopt.cellID{k}             = str2num((expcellids6{k}{1}));
 
 
   batchopt.loaddrive{k}        = xls_txt(i,loadrivecol);
 
    k = k+1;
end

