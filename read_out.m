%extract all cells in vivo parameters
load_depth=1;%load pial depth for all in vivo cells 
ODp=1;%Ocular dominance visual stimulation protocol 
SFTFp=1;
Sponp=1;
savefile=1;
%% 
%Excel sheet with all information
ExpXls            = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\L23PC_project_all_invivo.xlsx';%directory where excel batch file is located;change accordingly
%Structure in vitro
str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str';
%OD info
adata_dir         = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\OD-20200208T204106Z-001\';%data directory of raw data;change accordingly
%SFTF info
sftf_path         = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\SFTF\all';
sftf_path_resp    = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\SFTF';
%spontanues
spon_path         = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Spon_act';
%Pial depth info
roi_z             = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Pial_depth';

out_dir           = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\str_out';%data directory of extracted date;change accordingly
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
     L23_PC(adder).OD=OD;
     L23_PC(adder).SFTF=SFTF; 
     L23_PC(adder).spon=spon; 
     L23_PC(adder).pial_depth=pial_depth;
    adder=adder+1;
 end 
 
 %% 
 %% LOAD data structure, you need uipickfiles function
out_dir='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
 %% EXTRACT NECESSARY INFO from structure for iviv cells 
 for i=1:length(L23_PC)
     temp=L23_PC(i).ivivROI;
     oresp{i,:}=L23_PC(i).OD.oresp(temp);
     contra{i,:}=L23_PC(i).OD.contra(temp);
     ipsi{i,:}=L23_PC(i).OD.ipsi(temp);
     gOSI{i,:}=L23_PC(i).OD.gOSI(temp,:);
     gDSI{i,:}=L23_PC(i).OD.gDSI(temp,:);
     ODI{i,:}=L23_PC(i).OD.ODI(temp);
     allCa{i,:}=[sum(L23_PC(i).OD.deltapeaks_averagetrace(temp,1:8),2) sum(L23_PC(i).OD.deltapeaks_averagetrace(temp,9:16),2)];
     prefOri{i,:}=L23_PC(i).OD.prefOri(temp,:);
     prefDir{i,:}=L23_PC(i).OD.prefDir(temp,:);
     oppResp{i,:}=L23_PC(i).OD.oppResp(temp,:);
     sigma{i,:}=L23_PC(i).OD.sigma(temp,:);
     fit_tuning{i,:}=L23_PC(i).OD.fit_tuning(temp,:)';
     pia{i,:}=L23_PC(i).pial_depth(temp);
     temp=[];
 end
 %% 
 Dir=vertcat(prefDir{:});
 Ori=vertcat(prefOri{:});
 oppResp=vertcat(oppResp{:});
 sigma_tuning=vertcat(sigma{:});
 fit_resp=horzcat(fit_tuning{:});
 respo=horzcat(oresp{:});
 contra=horzcat(contra{:});
 ipsi=horzcat(ipsi{:});
 OSI=vertcat(gOSI{:});
 DSI=vertcat(gDSI{:});
 ODI=horzcat(ODI{:});
 Ca_sum=vertcat(allCa{:});
 pia_all=horzcat(pia{:});
 %% 
 %SFTF
   
      
 for i=1:length(L23_PC)
     temp=L23_PC(i).ivivROI;
  
     SFp{i,:}=L23_PC(i).SFTF.sf(temp);
     TFp{i,:}=L23_PC(i).SFTF.tf(temp);
     Ori_sftf{i,:}=L23_PC(i).SFTF.oripref(temp);
     Dir_sftf{i,:}=L23_PC(i).SFTF.dirpref(temp);
     Ca_peak{i,:}=L23_PC(i).SFTF.peakresp(temp);
     %for the different combination of SFTF
     Ori_sftf1{i,:}=L23_PC(i).SFTF.oripref1(temp);
     Dir_sftf1{i,:}=L23_PC(i).SFTF.dirpref1(temp);
     Ori_sftf2{i,:}=L23_PC(i).SFTF.oripref2(temp);
     Dir_sftf2{i,:}=L23_PC(i).SFTF.dirpref2(temp);
     Ori_sftf3{i,:}=L23_PC(i).SFTF.oripref3(temp);
     Dir_sftf3{i,:}=L23_PC(i).SFTF.dirpref3(temp);
      Ori_sftf4{i,:}=L23_PC(i).SFTF.oripref4(temp);
     Dir_sftf4{i,:}=L23_PC(i).SFTF.dirpref4(temp);
      Ori_sftf5{i,:}=L23_PC(i).SFTF.oripref5(temp);
     Dir_sftf5{i,:}=L23_PC(i).SFTF.dirpref5(temp);
      Ori_sftf6{i,:}=L23_PC(i).SFTF.oripref6(temp);
     Dir_sftf6{i,:}=L23_PC(i).SFTF.dirpref6(temp);
      Ori_sftf7{i,:}=L23_PC(i).SFTF.oripref7(temp);
     Dir_sftf7{i,:}=L23_PC(i).SFTF.dirpref7(temp);
      Ori_sftf8{i,:}=L23_PC(i).SFTF.oripref8(temp);
     Dir_sftf8{i,:}=L23_PC(i).SFTF.dirpref8(temp);
      Ori_sftf9{i,:}=L23_PC(i).SFTF.oripref9(temp);
     Dir_sftf9{i,:}=L23_PC(i).SFTF.dirpref9(temp);
     %OSI
     OSI_sftf{i,:}=L23_PC(i).SFTF.OSI(temp);
     DSI_sftf{i,:}=L23_PC(i).SFTF.DSI(temp);
     %OSI all combis
      OSI_sftf1{i,:}=L23_PC(i).SFTF.OSI1(temp);
     DSI_sftf1{i,:}=L23_PC(i).SFTF.DSI1(temp);
      OSI_sftf2{i,:}=L23_PC(i).SFTF.OSI2(temp);
     DSI_sftf2{i,:}=L23_PC(i).SFTF.DSI2(temp);
      OSI_sftf3{i,:}=L23_PC(i).SFTF.OSI3(temp);
     DSI_sftf3{i,:}=L23_PC(i).SFTF.DSI3(temp);
      OSI_sftf4{i,:}=L23_PC(i).SFTF.OSI4(temp);
     DSI_sftf4{i,:}=L23_PC(i).SFTF.DSI(temp);
      OSI_sftf5{i,:}=L23_PC(i).SFTF.OSI5(temp);
     DSI_sftf5{i,:}=L23_PC(i).SFTF.DSI5(temp);
      OSI_sftf6{i,:}=L23_PC(i).SFTF.OSI6(temp);
     DSI_sftf6{i,:}=L23_PC(i).SFTF.DSI6(temp);
      OSI_sftf7{i,:}=L23_PC(i).SFTF.OSI7(temp);
     DSI_sftf7{i,:}=L23_PC(i).SFTF.DSI7(temp);
      OSI_sftf8{i,:}=L23_PC(i).SFTF.OSI8(temp);
     DSI_sftf8{i,:}=L23_PC(i).SFTF.DSI8(temp);
      OSI_sftf9{i,:}=L23_PC(i).SFTF.OSI9(temp);
     DSI_sftf9{i,:}=L23_PC(i).SFTF.DSI9(temp);
     %responsive or not
     sftf_resp{i,:}=L23_PC(i).SFTF.ov_resp(temp);
     sftf_resp1{i,:}=L23_PC(i).SFTF.res(1,temp);
     sftf_resp2{i,:}=L23_PC(i).SFTF.res(2,temp);
     sftf_resp3{i,:}=L23_PC(i).SFTF.res(3,temp);
     sftf_resp4{i,:}=L23_PC(i).SFTF.res(4,temp);
     sftf_resp5{i,:}=L23_PC(i).SFTF.res(5,temp);
     sftf_resp6{i,:}=L23_PC(i).SFTF.res(6,temp);
     sftf_resp7{i,:}=L23_PC(i).SFTF.res(7,temp);
     sftf_resp8{i,:}=L23_PC(i).SFTF.res(8,temp);
     sftf_resp9{i,:}=L23_PC(i).SFTF.res(9,temp);
     %Peak Ca
     Ca_peak1{i,:}=L23_PC(i).SFTF.peakresp1(temp);
     Ca_peak2{i,:}=L23_PC(i).SFTF.peakresp2(temp);
     Ca_peak3{i,:}=L23_PC(i).SFTF.peakresp3(temp);
     Ca_peak4{i,:}=L23_PC(i).SFTF.peakresp4(temp);
     Ca_peak5{i,:}=L23_PC(i).SFTF.peakresp5(temp);
     Ca_peak6{i,:}=L23_PC(i).SFTF.peakresp6(temp);
     Ca_peak7{i,:}=L23_PC(i).SFTF.peakresp7(temp);
     Ca_peak8{i,:}=L23_PC(i).SFTF.peakresp8(temp);
     Ca_peak9{i,:}=L23_PC(i).SFTF.peakresp9(temp);
     
         for k=1:length(temp)
             if isempty(L23_PC(i).SFTF.Fit)==0
             if isempty(L23_PC(i).SFTF.Fit(temp(k)).PrefRsp)==0
%             Ori_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).PrefOri;
%             Dir_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).PrefDir;
            sigma_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).Sigma;
            oppResp_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).OppResp;
%             OSI_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).OSI;
%             DSI_sftf{i,k,:}=L23_PC(i).SFTF.Fit(temp(k)).DSI;
            %fit_sftf{i,k,:}=L23_PC(i).SFTF.Fit.FittedData(temp);
             else
%              Ori_sftf{i,k,:}=NaN;
%             Dir_sftf{i,k,:}=NaN;
            sigma_sftf{i,k,:}=NaN;;
            oppResp_sftf{i,k,:}=NaN;
%             OSI_sftf{i,k,:}=NaN;
%             DSI_sftf{i,k,:}=NaN;
         end
     else
%            Ori_sftf{i,k,:}=NaN;
%             Dir_sftf{i,k,:}=NaN;
            sigma_sftf{i,k,:}=NaN;;
            oppResp_sftf{i,k,:}=NaN;
%             OSI_sftf{i,k,:}=NaN;
%             DSI_sftf{i,k,:}=NaN;
             end
     
         end
         temp=[];
 end
 %% 
 Ori_sftf=horzcat(Ori_sftf{:});
 Dir_sftf=horzcat(Dir_sftf{:});
 
  Ori_sftf1=horzcat(Ori_sftf1{:});
 Dir_sftf1=horzcat(Dir_sftf1{:});
  Ori_sftf2=horzcat(Ori_sftf2{:});
 Dir_sftf2=horzcat(Dir_sftf2{:});
  Ori_sftf3=horzcat(Ori_sftf3{:});
 Dir_sftf3=horzcat(Dir_sftf3{:});
  Ori_sftf4=horzcat(Ori_sftf4{:});
 Dir_sftf4=horzcat(Dir_sftf4{:});
  Ori_sftf5=horzcat(Ori_sftf5{:});
 Dir_sftf5=horzcat(Dir_sftf5{:});
  Ori_sftf6=horzcat(Ori_sftf6{:});
 Dir_sftf6=horzcat(Dir_sftf6{:});
 Ori_sftf7=horzcat(Ori_sftf7{:});
 Dir_sftf7=horzcat(Dir_sftf7{:});
 Ori_sftf8=horzcat(Ori_sftf8{:});
 Dir_sftf8=horzcat(Dir_sftf8{:});
 Ori_sftf9=horzcat(Ori_sftf9{:});
 Dir_sftf9=horzcat(Dir_sftf9{:});
 
 oppResp_sftf=vertcat(oppResp_sftf{:});
 sigma_sftf=vertcat(sigma_sftf{:});
 sf=horzcat(SFp{:});
 tf=horzcat(TFp{:});
 sftf_resp=[sftf_resp{:}];
 
 OSI_sftf=horzcat(OSI_sftf{:});
 DSI_sftf=horzcat(DSI_sftf{:});
 Ca_peak_sf=horzcat(Ca_peak{:});
 
  OSI_sftf1=horzcat(OSI_sftf1{:});
 DSI_sftf1=horzcat(DSI_sftf1{:});
  OSI_sftf2=horzcat(OSI_sftf2{:});
 DSI_sftf2=horzcat(DSI_sftf2{:});
  OSI_sftf3=horzcat(OSI_sftf3{:});
 DSI_sftf3=horzcat(DSI_sftf3{:});
  OSI_sftf4=horzcat(OSI_sftf4{:});
 DSI_sftf4=horzcat(DSI_sftf4{:});
  OSI_sftf5=horzcat(OSI_sftf5{:});
 DSI_sftf5=horzcat(DSI_sftf5{:});
  OSI_sftf6=horzcat(OSI_sftf6{:});
 DSI_sftf6=horzcat(DSI_sftf6{:});
  OSI_sftf7=horzcat(OSI_sftf7{:});
 DSI_sftf7=horzcat(DSI_sftf7{:});
  OSI_sftf8=horzcat(OSI_sftf8{:});
 DSI_sftf8=horzcat(DSI_sftf8{:});
  OSI_sftf9=horzcat(OSI_sftf9{:});
 DSI_sftf9=horzcat(DSI_sftf9{:});
 
 Ca_peak_sf1=horzcat(Ca_peak1{:});
 Ca_peak_sf2=horzcat(Ca_peak2{:});
 Ca_peak_sf3=horzcat(Ca_peak3{:});
 Ca_peak_sf4=horzcat(Ca_peak4{:});
 Ca_peak_sf5=horzcat(Ca_peak5{:});
 Ca_peak_sf6=horzcat(Ca_peak6{:});
 Ca_peak_sf7=horzcat(Ca_peak7{:});
 Ca_peak_sf8=horzcat(Ca_peak8{:});
 Ca_peak_sf9=horzcat(Ca_peak9{:});
 
 sfre1=horzcat(sftf_resp1{:});
 sfre2=horzcat(sftf_resp2{:});
 sfre3=horzcat(sftf_resp3{:});
 sfre4=horzcat(sftf_resp4{:});
 sfre5=horzcat(sftf_resp5{:});
 sfre6=horzcat(sftf_resp6{:});
 sfre7=horzcat(sftf_resp7{:});
 sfre8=horzcat(sftf_resp8{:});
 sfre9=horzcat(sftf_resp9{:});

 
 %% Clear out nonsense from missing data: use OSI/DSI!
 
 
 idxn=find(isnan(OSI_sftf)==1);
 Ori_sftf(idxn)=NaN;
 Dir_sftf(idxn)=NaN;
 sigma_sftf(idxn)=NaN;
 Ca_peak_sf(idxn)=NaN;
 
 OSI_sftf1(find(sfre1==0))==NaN;
 OSI_sftf2(find(sfre2==0))==NaN;
 OSI_sftf3(find(sfre3==0))==NaN;
 OSI_sftf4(find(sfre4==0))==NaN;
 OSI_sftf5(find(sfre5==0))==NaN;
 OSI_sftf6(find(sfre6==0))==NaN;
 OSI_sftf7(find(sfre7==0))==NaN;
 OSI_sftf8(find(sfre8==0))==NaN;
 OSI_sftf9(find(sfre9==0))==NaN;
 
  DSI_sftf1(find(sfre1==0))==NaN;
 DSI_sftf2(find(sfre2==0))==NaN;
 DSI_sftf3(find(sfre3==0))==NaN;
DSI_sftf4(find(sfre4==0))==NaN;
 DSI_sftf5(find(sfre5==0))==NaN;
 DSI_sftf6(find(sfre6==0))==NaN;
 DSI_sftf7(find(sfre7==0))==NaN;
 DSI_sftf8(find(sfre8==0))==NaN;
 DSI_sftf9(find(sfre9==0))==NaN;
 
 idxn1=find(isnan(OSI_sftf1)==1);
 idxn2=find(isnan(OSI_sftf2)==1);
 idxn3=find(isnan(OSI_sftf3)==1);
 idxn4=find(isnan(OSI_sftf4)==1);
 idxn5=find(isnan(OSI_sftf5)==1);
 idxn6=find(isnan(OSI_sftf6)==1);
 idxn7=find(isnan(OSI_sftf7)==1);
 idxn8=find(isnan(OSI_sftf8)==1);
 idxn9=find(isnan(OSI_sftf9)==1);
  Ori_sftf1(idxn1)=NaN;
Dir_sftf1(idxn1)=NaN;
  Ca_peak_sf1(idxn1)=NaN;
  Ori_sftf2(idxn2)=NaN;
 Dir_sftf2(idxn2)=NaN;
 Ca_peak_sf2(idxn2)=NaN;
  Ori_sftf3(idxn3)=NaN;
 Dir_sftf3(idxn3)=NaN;
 Ca_peak_sf3(idxn3)=NaN;
  Ori_sftf4(idxn4)=NaN;
 Dir_sftf4(idxn4)=NaN;
 Ca_peak_sf4(idxn4)=NaN;
  Ori_sftf5(idxn5)=NaN;
 Dir_sftf5(idxn5)=NaN;
 Ca_peak_sf5(idxn5)=NaN;
  Ori_sftf6(idxn6)=NaN;
 Dir_sftf6(idxn6)=NaN;
  Ca_peak_sf6(idxn6)=NaN;
  Ori_sftf7(idxn7)=NaN;
 Dir_sftf7(idxn7)=NaN;
 Ca_peak_sf7(idxn7)=NaN;
  Ori_sftf8(idxn8)=NaN;
 Dir_sftf8(idxn8)=NaN;
  Ca_peak_sf8(idxn8)=NaN;
  Ori_sftf9(idxn9)=NaN;
 Dir_sftf9(idxn9)=NaN;
 Ca_peak_sf9(idxn9)=NaN;


 %% 
  for i=1:length(L23_PC)
      temp=L23_PC(i).ivivROI;
 pc{i,:}=L23_PC(i).spon.pci(temp);
 sa{i,:}=L23_PC(i).spon.sad(temp);
  end
  pci=vertcat(pc{:});
  sad=horzcat(sa{:});
  
  %% 
  str_invitro       = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% 

iviv_cells = [str(:).iviv]==1;
str_m=str;

aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
    if contra(i)==1 & ipsi(i)==0
        str_m(aiviv(i)).contra=1;
    elseif contra(i)==0 & ipsi(i)==1
   str_m(aiviv(i)).ipsi=1;
    elseif contra(i)==1 & ipsi(i)==1
     str_m(aiviv(i)).bino=1;   
    else contra(i)==0 & ipsi(i)==0
        str_m(aiviv(i)).unres=1;
       
    end
end
%% 

aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
   str_m(aiviv(i)).Ori=Ori(i,:);
   str_m(aiviv(i)).Dir=Dir(i,:);
   str_m(aiviv(i)).oppResp=oppResp(i,:);
   str_m(aiviv(i)).sigma=sigma_tuning(i,:);
   str_m(aiviv(i)).fit_resp=fit_resp(:,i);
   str_m(aiviv(i)).resp=respo(i);
   str_m(aiviv(i)).ODI=ODI(i);
   str_m(aiviv(i)).OSI=OSI(i,:);
   str_m(aiviv(i)).DSI=DSI(i,:);
   str_m(aiviv(i)).Ca_sum=Ca_sum(i,:);
   str_m(aiviv(i)).sftf_resp=sftf_resp(i);
   str_m(aiviv(i)).SF=sf(i);
   str_m(aiviv(i)).TF=tf(i);
   str_m(aiviv(i)).Ori_sftf=Ori_sftf(i);
   str_m(aiviv(i)).Dir_sftf=Dir_sftf(i); 
    str_m(aiviv(i)).Ca_peak_sftf=Ca_peak_sf(i); 
   str_m(aiviv(i)).Ori_sftf_all=[Ori_sftf1(i) Ori_sftf2(i) Ori_sftf3(i) Ori_sftf4(i) Ori_sftf5(i) Ori_sftf6(i) Ori_sftf7(i) Ori_sftf8(i) Ori_sftf9(i)];
   str_m(aiviv(i)).Dir_sftf_all=[Dir_sftf1(i) Dir_sftf2(i) Dir_sftf3(i) Dir_sftf4(i) Dir_sftf5(i) Dir_sftf6(i) Dir_sftf7(i) Dir_sftf8(i) Dir_sftf9(i)];
   str_m(aiviv(i)).reso_sftf_all=[sfre1(i) sfre2(i) sfre3(i) sfre4(i) sfre5(i) sfre6(i) sfre7(i) sfre8(i) sfre9(i)];
    str_m(aiviv(i)).Ca_sftf_all=[Ca_peak_sf1(i) Ca_peak_sf2(i) Ca_peak_sf3(i) Ca_peak_sf4(i) Ca_peak_sf5(i) Ca_peak_sf6(i) Ca_peak_sf7(i) Ca_peak_sf8(i) Ca_peak_sf9(i)];
   str_m(aiviv(i)).oppResp_sftf=oppResp_sftf(i); 
   str_m(aiviv(i)).sigma_sftf=sigma_sftf(i); 
   str_m(aiviv(i)).OSI_sftf=1-OSI_sftf(i); 
   str_m(aiviv(i)).DSI_sftf=1-DSI_sftf(i);
   str_m(aiviv(i)).OSI_sftf_all=[1-OSI_sftf1(i) 1-OSI_sftf2(i) 1-OSI_sftf3(i) 1-OSI_sftf4(i) 1-OSI_sftf5(i) 1-OSI_sftf6(i) 1-OSI_sftf7(i) 1-OSI_sftf8(i) 1-OSI_sftf9(i)];
   str_m(aiviv(i)).DSI_sftf_all=[1-DSI_sftf1(i) 1-DSI_sftf2(i) 1-DSI_sftf3(i) 1-DSI_sftf4(i) 1-DSI_sftf5(i) 1-DSI_sftf6(i) 1-DSI_sftf7(i) 1-DSI_sftf8(i) 1-DSI_sftf9(i)]; 
   str_m(aiviv(i)).pci=pci(i); 
   str_m(aiviv(i)).sad=sad(i); 
end


%%
aiviv=[];
aiviv=non_nan_idx;
for i=1:length(aiviv)
   str_m(aiviv(i)).ang_wmapi=ang_ai(i);
   str_m(aiviv(i)).ang_wmapil3=ang_ail3(i);
   str_m(aiviv(i)).ang_wmapil4=ang_ail4(i);
   
   str_m(aiviv(i)).ang_wmape=ang_ae(i);
   str_m(aiviv(i)).ang_wmapel3=ang_ael3(i);
   str_m(aiviv(i)).ang_wmapel4=ang_ael4(i);
end
%% 
aiviv=non_nan_idx;
for i=1:length(aiviv)
   str_m(aiviv(i)).Oripref=oripref(i);
 
end
%% 
aiviv=non_nan_idx;
for i=1:length(aiviv)
   str_m(aiviv(i)).PC=data_w_input(i,:);
 
end
%% 
iviv_cells = [str(:).iviv]==1;
str_m=str;
aiviv=find(iviv_cells==1);

for i=1:length(aiviv)
   str_m(aiviv(i)).pia_invivo=pia_all(i);

end
