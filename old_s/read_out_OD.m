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
     
      L23_PC(adder).mouse=[char(batchopt.mouse{i})];
     L23_PC(adder).exp_invitro=[char(batchopt.exp_invitro{i})];
     L23_PC(adder).ivivROI=batchopt.ivivROI{i};
     L23_PC(adder).OD=OD;
     adder=adder+1;
end
%% 
for i=1:31
L23_PC(i).OD.contra_tun=tt(i).OD.contra_tun
L23_PC(i).OD.ipsi_tun=tt(i).OD.ipsi_tun
L23_PC(i).OD.fit_error=tt(i).OD.fit_error
L23_PC(i).OD.fit_r2=tt(i).OD.fit_r2
end
