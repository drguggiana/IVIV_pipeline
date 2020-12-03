%extract all cells in vivo parameters
load_depth=1;%load pial depth for all in vivo cells 
ODp=1;%Ocular dominance visual stimulation protocol 
SFTFp=1;
Sponp=1;
savefile=1;
%% 
%Excel sheet with all information
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\L23PC_project_all_invivo.xlsx';%directory where excel batch file is located;change accordingly
%Structure in vitro
str_invitro       = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\old_structures';
%OD info
adata_dir         = 'D:\Postdoc_Margrie\Projects\L23\data_invivo\OD\';%data directory of raw data;change accordingly
%SFTF info
sftf_path         = 'D:\Postdoc_Margrie\Projects\L23\data_invivo\SFTF\all';
sftf_path_resp    = 'D:\Postdoc_Margrie\Projects\L23\data_invivo\SFTF';
%spontanues
spon_path         = 'D:\Postdoc_Margrie\Projects\L23\data_invivo\Spon_act';
%Pial depth info
roi_z             = 'D:\Postdoc_Margrie\Projects\L23\data_invivo\pial_depth_output';

out_dir           = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo';%data directory of extracted date;change accordingly
%%   Read out excel info
batchopt          = parseExperimentsXls_L23(ExpXls);%calls the nested function parseExperimentsXls_dLGN and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
%% 
 adder=1;%counting variable
 for i=1:nummice%for loop over experiments across days
     %OD
     pathName=adata_dir;
     fileOD=dir([char(pathName) '\exp' mat2str(batchopt.binoexp_ids{i}) '\ORIanalyis_z8_meds_trace_mean_zeroneg_fit_NP_detrend_globalF0.mat']);
     cd(char([char(pathName) '\exp'   mat2str(batchopt.binoexp_ids{i})]));
     load(fileOD.name);
     [OD] = extractOD_mod(peaks,Fit);
     peaks=[];
     Fit=[];
     %sftf resp
     fileSFTF=dir([char(sftf_path_resp) '\block_exp' mat2str(batchopt.sftfexp_ids{i}(1)) '_to_'  mat2str(batchopt.sftfexp_ids{i}(end)) '_allSFTF_consolidated.mat']);
     cd(fileSFTF.folder)
     load(fileSFTF.name);
     [SFTF] = extractSFTF(peaks);
    % normr_2(permute(peaks.average_trace_delta_integ,[1 3 4 2])).*255
     %sftf data
     SF=[0.02 0.08 0.16];
     TF=[1 2 4];
     oris=[0;45;90;135;180;225;270;315];
     cd(sftf_path);
     load('180601_1649_sftfData_FIXED.mat');
     temp=sftf_mat(i).data;
      
for k=1:length(temp)
     if ~isnan(temp(k,1));
       temp1=squeeze(nanmean(max(temp(k,:,:,:),[],4),1));
        [M,I] = max(temp1(:));
     [I_SF, I_TF] = ind2sub(size(temp1),I);
     tf_c(k)=TF(I_TF);
     sf_c(k)=SF(I_SF);
     %orie fit at best SF TF
     oristft=squeeze(temp(k,I_SF,I_TF,1:8))';
    % figure;plot(oris,oristft);
     [Fit(k).FittedData, Fit(k).BaselineRsp, Fit(k).PrefRsp, Fit(k).PrefDir, Fit(k).Sigma, Fit(k).OppResp ,...
      Fit(k).Error, Fit(k).R2]=FIT_Carandini(oristft);
      Fit(k).orifit= TT_FoldedTuningCurve(Fit(k).FittedData);
      Fit(k).PrefOri = [mod([Fit(k).PrefDir]+180-1, 180)+1]';
      OSI(k)=TT_CircularVariance_ORI(oristft);
      DSI(k)=TT_CircularVariance(oristft);
      oripref(k)=Fit(k).PrefOri;
      dirpref(k)=Fit(k).PrefDir;
      peakresp(k)=M;
      % for individual SF TF combis
      %#1
      oritemp=squeeze(temp(k,1,1,1:8))';
      [Fit1(k).FittedData, Fit1(k).BaselineRsp, Fit1(k).PrefRsp, Fit1(k).PrefDir, Fit1(k).Sigma, Fit1(k).OppResp ,...
      Fit1(k).Error, Fit1(k).R2]=FIT_Carandini(oritemp);
      Fit1(k).orifit= TT_FoldedTuningCurve(Fit1(k).FittedData);
      Fit1(k).PrefOri = [mod([Fit1(k).PrefDir]+180-1, 180)+1]';
      OSI1(k)=TT_CircularVariance_ORI(oritemp);
      DSI1(k)=TT_CircularVariance(oritemp);
      oripref1(k)=Fit1(k).PrefOri;
      dirpref1(k)=Fit1(k).PrefDir;
      peakresp1(k)=max(oritemp);
      oritemp=[];
      %#2
      oritemp=squeeze(temp(k,2,1,1:8))';
      [Fit2(k).FittedData, Fit2(k).BaselineRsp, Fit2(k).PrefRsp, Fit2(k).PrefDir, Fit2(k).Sigma, Fit2(k).OppResp ,...
      Fit2(k).Error, Fit2(k).R2]=FIT_Carandini(oritemp);
      Fit2(k).orifit= TT_FoldedTuningCurve(Fit2(k).FittedData);
      Fit2(k).PrefOri = [mod([Fit2(k).PrefDir]+180-1, 180)+1]';
      OSI2(k)=TT_CircularVariance_ORI(oritemp);
      DSI2(k)=TT_CircularVariance(oritemp);
      oripref2(k)=Fit2(k).PrefOri;
      dirpref2(k)=Fit2(k).PrefDir;
      peakresp2(k)=max(oritemp);
      oritemp=[];
      %#3
      oritemp=squeeze(temp(k,3,1,1:8))';
      [Fit3(k).FittedData, Fit3(k).BaselineRsp, Fit3(k).PrefRsp, Fit3(k).PrefDir, Fit3(k).Sigma, Fit3(k).OppResp ,...
      Fit3(k).Error, Fit3(k).R2]=FIT_Carandini(oritemp);
      Fit3(k).orifit= TT_FoldedTuningCurve(Fit3(k).FittedData);
      Fit3(k).PrefOri = [mod([Fit3(k).PrefDir]+180-1, 180)+1]';
      OSI3(k)=TT_CircularVariance_ORI(oritemp);
      DSI3(k)=TT_CircularVariance(oritemp);
      oripref3(k)=Fit3(k).PrefOri;
      dirpref3(k)=Fit3(k).PrefDir;
      peakresp3(k)=max(oritemp);
      oritemp=[];
      %#4
      oritemp=squeeze(temp(k,1,2,1:8))';
      [Fit4(k).FittedData, Fit4(k).BaselineRsp, Fit4(k).PrefRsp, Fit4(k).PrefDir, Fit4(k).Sigma, Fit4(k).OppResp ,...
      Fit4(k).Error, Fit4(k).R2]=FIT_Carandini(oritemp);
      Fit4(k).orifit= TT_FoldedTuningCurve(Fit4(k).FittedData);
      Fit4(k).PrefOri = [mod([Fit4(k).PrefDir]+180-1, 180)+1]';
      OSI4(k)=TT_CircularVariance_ORI(oritemp);
      DSI4(k)=TT_CircularVariance(oritemp);
      oripref4(k)=Fit4(k).PrefOri;
      dirpref4(k)=Fit4(k).PrefDir;
      peakresp4(k)=max(oritemp);
      oritemp=[];
        %#5
      oritemp=squeeze(temp(k,2,2,1:8))';
      [Fit5(k).FittedData, Fit5(k).BaselineRsp, Fit5(k).PrefRsp, Fit5(k).PrefDir, Fit5(k).Sigma, Fit5(k).OppResp ,...
      Fit5(k).Error, Fit5(k).R2]=FIT_Carandini(oritemp);
      Fit5(k).orifit= TT_FoldedTuningCurve(Fit5(k).FittedData);
      Fit5(k).PrefOri = [mod([Fit5(k).PrefDir]+180-1, 180)+1]';
      OSI5(k)=TT_CircularVariance_ORI(oritemp);
      DSI5(k)=TT_CircularVariance(oritemp);
      oripref5(k)=Fit5(k).PrefOri;
      dirpref5(k)=Fit5(k).PrefDir;
      peakresp5(k)=max(oritemp);
      oritemp=[];
       %#6
      oritemp=squeeze(temp(k,3,2,1:8))';
      [Fit6(k).FittedData, Fit6(k).BaselineRsp, Fit6(k).PrefRsp, Fit6(k).PrefDir, Fit6(k).Sigma, Fit6(k).OppResp ,...
      Fit6(k).Error, Fit6(k).R2]=FIT_Carandini(oritemp);
      Fit6(k).orifit= TT_FoldedTuningCurve(Fit6(k).FittedData);
      Fit6(k).PrefOri = [mod([Fit6(k).PrefDir]+180-1, 180)+1]';
      OSI6(k)=TT_CircularVariance_ORI(oritemp);
      DSI6(k)=TT_CircularVariance(oritemp);
      oripref6(k)=Fit6(k).PrefOri;
      dirpref6(k)=Fit6(k).PrefDir;
      peakresp6(k)=max(oritemp);
      oritemp=[];
        %#7
      oritemp=squeeze(temp(k,1,3,1:8))';
      [Fit7(k).FittedData, Fit7(k).BaselineRsp, Fit7(k).PrefRsp, Fit7(k).PrefDir, Fit7(k).Sigma, Fit7(k).OppResp ,...
      Fit7(k).Error, Fit7(k).R2]=FIT_Carandini(oritemp);
      Fit7(k).orifit= TT_FoldedTuningCurve(Fit7(k).FittedData);
      Fit7(k).PrefOri = [mod([Fit7(k).PrefDir]+180-1, 180)+1]';
      OSI7(k)=TT_CircularVariance_ORI(oritemp);
      DSI7(k)=TT_CircularVariance(oritemp);
      oripref7(k)=Fit7(k).PrefOri;
      dirpref7(k)=Fit7(k).PrefDir;
      peakresp7(k)=max(oritemp);
      oritemp=[];
        %#8
      oritemp=squeeze(temp(k,2,3,1:8))';
      [Fit8(k).FittedData, Fit8(k).BaselineRsp, Fit8(k).PrefRsp, Fit8(k).PrefDir, Fit8(k).Sigma, Fit8(k).OppResp ,...
      Fit8(k).Error, Fit8(k).R2]=FIT_Carandini(oritemp);
      Fit8(k).orifit= TT_FoldedTuningCurve(Fit8(k).FittedData);
      Fit8(k).PrefOri = [mod([Fit8(k).PrefDir]+180-1, 180)+1]';
      OSI8(k)=TT_CircularVariance_ORI(oritemp);
      DSI8(k)=TT_CircularVariance(oritemp);
      oripref8(k)=Fit8(k).PrefOri;
      dirpref8(k)=Fit8(k).PrefDir;
      peakresp8(k)=max(oritemp);
      oritemp=[];
       %#9
      oritemp=squeeze(temp(k,3,3,1:8))';
      [Fit9(k).FittedData, Fit9(k).BaselineRsp, Fit9(k).PrefRsp, Fit9(k).PrefDir, Fit9(k).Sigma, Fit9(k).OppResp ,...
      Fit9(k).Error, Fit9(k).R2]=FIT_Carandini(oritemp);
      Fit9(k).orifit= TT_FoldedTuningCurve(Fit9(k).FittedData);
      Fit9(k).PrefOri = [mod([Fit9(k).PrefDir]+180-1, 180)+1]';
      OSI9(k)=TT_CircularVariance_ORI(oritemp);
      DSI9(k)=TT_CircularVariance(oritemp);
      oripref9(k)=Fit9(k).PrefOri;
      dirpref9(k)=Fit9(k).PrefDir;
      peakresp9(k)=max(oritemp);
      oritemp=[];
     else
        tf_c(k)=NaN;
        sf_c(k)=NaN;
        oristft=NaN;
        OSI(k)=NaN;OSI1(k)=NaN;OSI2(k)=NaN;OSI3(k)=NaN;OSI4(k)=NaN;OSI5(k)=NaN;OSI6(k)=NaN;OSI7(k)=NaN;OSI8(k)=NaN;OSI9(k)=NaN;
        DSI(k)=NaN;DSI1(k)=NaN;DSI2(k)=NaN;DSI3(k)=NaN;DSI4(k)=NaN;DSI5(k)=NaN;DSI6(k)=NaN;DSI7(k)=NaN;DSI8(k)=NaN;DSI9(k)=NaN;
        oripref(k)=NaN;oripref1(k)=NaN;oripref2(k)=NaN;oripref3(k)=NaN;oripref4(k)=NaN;oripref5(k)=NaN;oripref6(k)=NaN;oripref7(k)=NaN;oripref8(k)=NaN;oripref9(k)=NaN;
        dirpref(k)=NaN;dirpref1(k)=NaN;dirpref2(k)=NaN;dirpref3(k)=NaN;dirpref4(k)=NaN;dirpref5(k)=NaN;dirpref6(k)=NaN;dirpref7(k)=NaN;dirpref8(k)=NaN;dirpref9(k)=NaN;
        peakresp1(k)=NaN;peakresp2(k)=NaN;peakresp3(k)=NaN;peakresp4(k)=NaN;peakresp5(k)=NaN;peakresp6(k)=NaN;peakresp7(k)=NaN;peakresp8(k)=NaN;peakresp9(k)=NaN;
        peakresp(k)=NaN;
     end
     
end

SFTF.expname= sftf_mat(i).expName;
SFTF.data=  sftf_mat(i).data;
SFTF.sf=sf_c;
SFTF.tf=tf_c;
SFTF.Fit=Fit;SFTF.Fit1=Fit1;SFTF.Fit2=Fit2;SFTF.Fit3=Fit3;SFTF.Fit4=Fit4;SFTF.Fit5=Fit5;SFTF.Fit6=Fit6;SFTF.Fit7=Fit7;SFTF.Fit8=Fit8;SFTF.Fit9=Fit9;
SFTF.OSI=OSI;SFTF.OSI1=OSI1;SFTF.OSI2=OSI2;SFTF.OSI3=OSI3;SFTF.OSI4=OSI4;SFTF.OSI5=OSI5;SFTF.OSI6=OSI6;SFTF.OSI7=OSI7;SFTF.OSI8=OSI8;SFTF.OSI9=OSI9;
SFTF.DSI=DSI;SFTF.DSI1=DSI1;SFTF.DSI2=DSI2;SFTF.DSI3=DSI3;SFTF.DSI4=DSI4;SFTF.DSI5=DSI5;SFTF.DSI6=DSI6;SFTF.DSI7=DSI7;SFTF.DSI8=DSI8;SFTF.DSI9=DSI9;
SFTF.oripref=oripref;SFTF.oripref1=oripref1;SFTF.oripref2=oripref2;SFTF.oripref3=oripref3;SFTF.oripref4=oripref4;SFTF.oripref5=oripref5;SFTF.oripref6=oripref6;
SFTF.oripref7=oripref7;SFTF.oripref8=oripref8;SFTF.oripref9=oripref9;
SFTF.dirpref=dirpref;SFTF.dirpref1=dirpref1;SFTF.dirpref2=dirpref2;SFTF.dirpref3=dirpref3;SFTF.dirpref4=dirpref4;SFTF.dirpref5=dirpref5;SFTF.dirpref6=dirpref6;
SFTF.dirpref7=dirpref7;SFTF.dirpref8=dirpref8;SFTF.dirpref9=dirpref9;
SFTF.peakresp=peakresp;SFTF.peakresp1=peakresp1;SFTF.peakresp2=peakresp2;SFTF.peakresp3=peakresp3;SFTF.peakresp4=peakresp4;SFTF.peakresp5=peakresp5;SFTF.peakresp6=peakresp6;
SFTF.peakresp7=peakresp7;SFTF.peakresp8=peakresp8;SFTF.peakresp9=peakresp9;
temp=[];
sf_c=[];
tf_c=[];
oristft=[];
Fit=[];Fit1=[];Fit2=[];Fit3=[];Fit4=[];Fit5=[];Fit6=[];Fit7=[];Fit8=[];Fit9=[];
OSI=[];OSI1=[];OSI2=[];OSI3=[];OSI4=[];OSI5=[];OSI6=[];OSI7=[];OSI8=[];OSI9=[];
DSI=[];DSI1=[];DSI2=[];DSI3=[];DSI4=[];DSI5=[];DSI6=[];DSI7=[];DSI8=[];DSI9=[];
oripref=[];oripref1=[];oripref2=[];oripref3=[];oripref4=[];oripref5=[];oripref6=[];oripref7=[];oripref8=[];oripref9=[];
dirpref=[];dirpref1=[];dirpref2=[];dirpref3=[];dirpref4=[];dirpref5=[];dirpref6=[];dirpref7=[];dirpref8=[];dirpref9=[];
peakresp=[];peakresp1=[];peakresp2=[];peakresp3=[];peakresp4=[];peakresp5=[];peakresp6=[];peakresp7=[];peakresp8=[];peakresp9=[];
%Load and extract spontaneous
cd(spon_path);
load('190606_1359_spontAct.mat');
idx=[1:30 33];
spon.pci=popcop_cell{idx(i)};
spon.sad=event_store{idx(i)};
%Load and extract pial depth 
   pathName_z=fullfile(roi_z , filesep);
   file_depth=dir([char(pathName_z) char(batchopt.mouse{i}) '_roiDepth.mat']);
   cd(char([char(pathName_z)]));
   load(file_depth.name);
   pial_depth=[rois(:).depth];  
     
     L23_PC(adder).mouse=[char(batchopt.mouse{i})];
     L23_PC(adder).exp_invitro=[char(batchopt.exp_invitro{i})];
     L23_PC(adder).ivivROI=batchopt.ivivROI{i};
     %L23_PC(adder).exp_invitro=[char(batchopt.invitro_rec{i})];
     L23_PC(adder).OD=OD;
     L23_PC(adder).SFTF=SFTF; 
     L23_PC(adder).spon=spon; 
     L23_PC(adder).pial_depth=pial_depth;
    adder=adder+1;
 end 
 
 %% Add other information 
 
