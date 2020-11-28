%% 
%% 
%1. You need the in vivo data structure called L23_PC which contains all in vivo cells
%and their extracted values from the OD data, SFTF data and spon activity
%data 

%%2. You need the invitro structure 

 %% LOAD in vivo data structure L23PC data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
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
     c_tun{i,:}=L23_PC(i).OD.contra_tun(temp);
     i_tun{i,:}=L23_PC(i).OD.ipsi_tun(temp);
     
     allCa{i,:}=[sum(L23_PC(i).OD.deltapeaks_averagetrace(temp,1:8),2) sum(L23_PC(i).OD.deltapeaks_averagetrace(temp,9:16),2)];
    peakCaOD{i,:}=[nanmax(L23_PC(i).OD.deltapeaks_averagetrace(temp,1:8),[],2) nanmax(L23_PC(i).OD.deltapeaks_averagetrace(temp,9:16),[],2)];
     prefOri{i,:}=L23_PC(i).OD.prefOri(temp,:);
     prefDir{i,:}=L23_PC(i).OD.prefDir(temp,:);
     oppResp{i,:}=L23_PC(i).OD.oppResp(temp,:);
     sigma{i,:}=L23_PC(i).OD.sigma(temp,:);
     erro{i,:}=L23_PC(i).OD.fit_error(temp,:);
     r2f{i,:}=L23_PC(i).OD.fit_r2(temp,:);
     fit_tuning{i,:}=L23_PC(i).OD.fit_tuning(temp,:)';
     pia{i,:}=L23_PC(i).pial_depth(temp);
     temp=[];
 end
 %% read out PSTH
 for i=1:length(L23_PC)
     temp=L23_PC(i).ivivROI;
     PSTH{i,:}=L23_PC(i).OD.deltapeaks_averagetrace(temp,1:16);
 end
 PSTH_all=vertcat(PSTH{:});
 
 %% concatenate cells OD protocol
 Dir=vertcat(prefDir{:});
 Ori=vertcat(prefOri{:});
 oppResp=vertcat(oppResp{:});
 sigma_tuning=vertcat(sigma{:});
 error_fit=vertcat(erro{:});
 r2_fit=vertcat(r2f{:});
 fit_resp=horzcat(fit_tuning{:});
 respo=horzcat(oresp{:});
 contra=horzcat(contra{:});
 ipsi=horzcat(ipsi{:});
 contra_tun=horzcat(c_tun{:});
 ipsi_tun=horzcat(i_tun{:});
 OSI=vertcat(gOSI{:});
 DSI=vertcat(gDSI{:});
 ODI=horzcat(ODI{:});
 Ca_sum=vertcat(allCa{:});
 Ca_peakOD=vertcat(peakCaOD{:});
 pia_all=horzcat(pia{:});
 %% concatenate cells SFTF protocol     
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
 %% Concatnate SFTF cells
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
 %% Spontaneuos activity: read out and concatenate
  for i=1:length(L23_PC)
      temp=L23_PC(i).ivivROI;
 pc{i,:}=L23_PC(i).spon.pci(temp);
 sa{i,:}=L23_PC(i).spon.sad(temp);
  end
  pci=vertcat(pc{:});
  sad=horzcat(sa{:});
 
 %% ASSEMBLE IVIV STRUCTURE
 %% 
 %% 
 
 
  %% Load in vitro strcuture
  str_invitro       = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\Stage1_invitro_only_structure';
folder_list = uipickfiles('FilterSpec',str_invitro);
load(char(folder_list));
%% 
  str_invitro2       = 'D:\Postdoc_Margrie\Projects\L23\structure';
folder_list = uipickfiles('FilterSpec',str_invitro2);
load(char(folder_list));


%% Read out iviv cells
iviv_cells = [str(:).iviv]==1;

%% 
%invitro_struct([44 48 68 71 91 94 109 121 156 157])=[]
%% 
% for i=1:147
% invitro_struct(i).morph=str(i).morph;
% invitro_struct(i).morphtraces=str(i).morphtraces;
% invitro_struct(i).iviv=str(i).iviv;
% invitro_struct(i).morphoMap_apical_aligned=str(i).morphoMap_apical_aligned;
% invitro_struct(i).morphoMap_basal_aligned=str(i).morphoMap_basal_aligned;
% end
%% 
str_m=invitro_struct;
%% 
for i=1:157
    str_m(i).iviv=iviv_cells(i);
end

%% Add contra, bino, ipsi, unresposnive for OD protocol
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
%% Replace empty fields with NaNs
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.ipsi),str_m))
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).ipsi=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.contra),str_m))
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).contra=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.bino),str_m))
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).bino=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.unres),str_m))
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).unres=NaN; 
end
%% ADD invivo features to strcutre, for OD protocol: first column contra, second column ipsi
aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
    %OD protocol
   str_m(aiviv(i)).Ori=Ori(i,:);
   str_m(aiviv(i)).Dir=Dir(i,:);
   %str_m(aiviv(i)).oppResp=oppResp(i,:);
   str_m(aiviv(i)).sigma=sigma_tuning(i,:);
   str_m(aiviv(i)).fit_resp=fit_resp(:,i);
   str_m(aiviv(i)).resp=respo(i);
   str_m(aiviv(i)).ODI=ODI(i);
   str_m(aiviv(i)).OSI=OSI(i,:);
   str_m(aiviv(i)).DSI=DSI(i,:);
   str_m(aiviv(i)).Ca_sum=Ca_sum(i,:);
   str_m(aiviv(i)).Ca_peak=Ca_peakOD(i,:);
   str_m(aiviv(i)).TuningCurve=PSTH_all(i,:);
   %SFTF protocol
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
   %str_m(aiviv(i)).oppResp_sftf=oppResp_sftf(i); 
   str_m(aiviv(i)).sigma_sftf=sigma_sftf(i); 
   str_m(aiviv(i)).OSI_sftf=1-OSI_sftf(i); 
   str_m(aiviv(i)).DSI_sftf=1-DSI_sftf(i);
   str_m(aiviv(i)).OSI_sftf_all=[1-OSI_sftf1(i) 1-OSI_sftf2(i) 1-OSI_sftf3(i) 1-OSI_sftf4(i) 1-OSI_sftf5(i) 1-OSI_sftf6(i) 1-OSI_sftf7(i) 1-OSI_sftf8(i) 1-OSI_sftf9(i)];
   str_m(aiviv(i)).DSI_sftf_all=[1-DSI_sftf1(i) 1-DSI_sftf2(i) 1-DSI_sftf3(i) 1-DSI_sftf4(i) 1-DSI_sftf5(i) 1-DSI_sftf6(i) 1-DSI_sftf7(i) 1-DSI_sftf8(i) 1-DSI_sftf9(i)]; 
   %spon activity protocol
   str_m(aiviv(i)).pci=pci(i); 
   str_m(aiviv(i)).sad=sad(i); 
   %pia from in vivo
   str_m(aiviv(i)).pia_invivo=pia_all(i);
end
%% 
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.Ori),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).Ori=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.Dir),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).Dir=NaN; 
end

emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.sigma),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).sigma=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.fit_resp),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).fit_resp=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.resp),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).resp=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.OSI),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).OSI=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.DSI),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).DSI=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.ODI),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).ODI=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.Ca_peak),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).Ca_peak=NaN; 
end

emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.Ca_sum),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).Ca_sum=NaN; 
end
emptyIndex = [];
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.TuningCurve),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).TuningCurve=NaN; 
end

emptyIndex = find(arrayfun(@(str_m) isempty(str_m.SF),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).SF=NaN; 
end

emptyIndex = find(arrayfun(@(str_m) isempty(str_m.TF),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).TF=NaN; 
end

emptyIndex = find(arrayfun(@(str_m) isempty(str_m.sftf_resp),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).sftf_resp=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.pci),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).pci=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.sad),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).sad=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.pia_invivo),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).pia_invivo=NaN; 
end
%% 
aiviv=find(iviv_cells==1);
for i=1:length(aiviv)
    str_m(aiviv(i)).c_tun=contra_tun(i);
    str_m(aiviv(i)).i_tun=ipsi_tun(i);
    str_m(aiviv(i)).error_fit=error_fit(i,:);
    str_m(aiviv(i)).r2_fit=r2_fit(i,:);
    
end
%% 
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.c_tun),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).c_tun=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.i_tun),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).i_tun=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.error_fit),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).error_fit=NaN; 
end
emptyIndex = find(arrayfun(@(str_m) isempty(str_m.r2_fit),str_m));
for i=1:length(emptyIndex)
   str_m(emptyIndex(i)).r2_fit=NaN; 
end
%% 
str_m([44 48 68 71 91 94 109 121 156 157])=[]
%% 

  str_invitro2       = 'D:\Postdoc_Margrie\Projects\L23\structure';
folder_list = uipickfiles('FilterSpec',str_invitro2);
load(char(folder_list));
%% 
 for i=1:147
 str_m(i).morph=str(i).morph;
 str_m(i).morphtraces=str(i).morphtraces;
 
 str_m(i).morphoMap_apical_aligned=str(i).morphoMap_apical_aligned;
str_m(i).morphoMap_basal_aligned=str(i).morphoMap_basal_aligned;
 end
%% 

str=str_m;
%% END FOR NOW

%% Extract eye specific info and add to strcuture 
%Here strcuture with 147 is necessary

%% 
%% Define array names and get the arrays without NaNs etc.
field_list = {'excMap','inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
iviv_cells = [str(:).iviv]==1;
morpho_cells = ~cellfun(@isempty, {str.morph});
% cell_idx = find(iviv_cells&morpho_cells);
% cell_idx = find(iviv_cells);
% cell_idx = find(morpho_cells);
cell_idx = 1:length(str);
% get the number of cells to include
cell_num = length(cell_idx);
% allocate memory for the individual cells
cell_cell = cell(cell_num,1);
% for all the cells
for cells = 1:cell_num
    % allocate a small cell to hold each field
    field_cell = cell(field_number,1);
    field_cell_raw = cell(field_number,1);
    
    % for all the fields
    for fields = 1:field_number
        field_cell{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'excMap')
            field_cell{fields} = field_cell{fields}./min(field_cell{fields});
        elseif contains(field_list{fields}, 'inhMap')
            field_cell{fields} = field_cell{fields}./max(field_cell{fields});
        end       
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});    
    
    for fields = 1:field_number
        field_cell_raw{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'excMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        elseif contains(field_list{fields}, 'inhMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        end       
    end
    
    cell_cell_raw{cells} = vertcat(field_cell_raw{:});
end
% concatenate the results
cell_cell = cat(2,cell_cell{:})';
cell_cell_raw = cat(2,cell_cell_raw{:})';
% remove cells with NaNs
non_nan_cells = sum(isnan(cell_cell),2)==0;
nan_vector = find(non_nan_cells>0);
cell_cell = cell_cell(non_nan_cells,:);
cell_cell_raw = cell_cell_raw(non_nan_cells,:);
% copy the cell to have the original maps later
original_maps = cell_cell;
original_maps_raw = cell_cell_raw;
% redefine cell number based on the rows that didn't contain NaN
cell_num = size(cell_cell,1);
%Pia vector for ex and inh maps 
pia_input=original_maps(:,end);
incl_idx=1;
%% Read out eye specific info

[od_out_iviv spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(str_m,1:147)
%% Add OD eye specific info to strcuture
nan_vector=1:147;
for i=1:length(nan_vector)
    str(nan_vector(i)).OSIpref=od_out_iviv(i,1);
    str(nan_vector(i)).DSIpref=od_out_iviv(i,2);
    str(nan_vector(i)).ODIpref=od_out_iviv(i,3);
    str(nan_vector(i)).ORIpref=od_out_iviv(i,4);
    str(nan_vector(i)).DIRpref=od_out_iviv(i,5);
    str(nan_vector(i)).Casumpref=od_out_iviv(i,6);
    str(nan_vector(i)).Sigmapref=od_out_iviv(i,7);
    str(nan_vector(i)).Capeakpref=od_out_iviv(i,8);
    str(nan_vector(i)).error_pref=od_out_iviv(i,10);
    str(nan_vector(i)).r2_pref=od_out_iviv(i,11);
    str(nan_vector(i)).tun_pref=od_out_iviv(i,12);
end
%% 
%% ADD preferred Tuning Curve
non_nan_idx=nan_vector(1:end);
for i=1:length(non_nan_idx);
%    idx_amp_inf(i,:)=str(non_nan_idx(i)).hori_peak_pl;
if ~isnan(str(non_nan_idx(i)).Ori)==1 & str(non_nan_idx(i)).resp==1
    
    if str(non_nan_idx(i)).contra==1
     tc(i,:)=str(non_nan_idx(i)).TuningCurve(1:8); 
     
    elseif str(non_nan_idx(i)).ipsi==1;
      tc(i,:)=str(non_nan_idx(i)).TuningCurve(9:end); 
    
      
    elseif str(non_nan_idx(i)).bino==1;
         if str(non_nan_idx(i)).ODI>=0;
        tc(i,:)=str(non_nan_idx(i)).TuningCurve(1:8); 
      
       
         else str(non_nan_idx(i)).ODI<0;
            tc(i,:)=str(non_nan_idx(i)).TuningCurve(9:end); 
      
       
    end
 else str(non_nan_idx(i)).unres==1;
        tc(i,:)=NaN*ones(1,8);
      
     
 end

else
tc(i,:)=NaN*ones(1,8);

end
end
%% 
for i=1:length(nan_vector)
    str(nan_vector(i)).TCpref=tc(i,:);
 
end

















%% Define array names and get the arrays without NaNs etc.
field_list = {'subpixel_excMap','subpixel_inhMap','pialD'};
% get the number of fields
field_number = length(field_list);
% get a vector with the cells to use
iviv_cells = [str(:).iviv]==1;
morpho_cells = ~cellfun(@isempty, {str.morph});
% cell_idx = find(iviv_cells&morpho_cells);
% cell_idx = find(iviv_cells);
% cell_idx = find(morpho_cells);
cell_idx = 1:length(str);
% get the number of cells to include
cell_num = length(cell_idx);
% allocate memory for the individual cells
cell_cell = cell(cell_num,1);
% for all the cells
for cells = 1:cell_num
    % allocate a small cell to hold each field
    field_cell = cell(field_number,1);
    field_cell_raw = cell(field_number,1);
    
    % for all the fields
    for fields = 1:field_number
        field_cell{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'subpixel_excMap')
            field_cell{fields} = field_cell{fields}./min(field_cell{fields});
        elseif contains(field_list{fields}, 'subpixel_inhMap')
            field_cell{fields} = field_cell{fields}./max(field_cell{fields});
        end       
    end
    % save the fields in the target cell
    cell_cell{cells} = vertcat(field_cell{:});    
    
    for fields = 1:field_number
        field_cell_raw{fields} = str(cell_idx(cells)).(field_list{fields})(:);
        if contains(field_list{fields}, 'subpixel_excMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        elseif contains(field_list{fields}, 'subpixel_inhMap')
            field_cell_raw{fields} = field_cell_raw{fields};
        end       
    end
    
    cell_cell_raw{cells} = vertcat(field_cell_raw{:});
end
% concatenate the results
cell_cell = cat(2,cell_cell{:})';
cell_cell_raw = cat(2,cell_cell_raw{:})';
% remove cells with NaNs
non_nan_cells = sum(isnan(cell_cell),2)==0;
nan_vector = find(non_nan_cells>0);
cell_cell = cell_cell(non_nan_cells,:);
cell_cell_raw = cell_cell_raw(non_nan_cells,:);
% copy the cell to have the original maps later
original_maps = cell_cell;
original_maps_raw = cell_cell_raw;
% redefine cell number based on the rows that didn't contain NaN
cell_num = size(cell_cell,1);
%Pia vector for ex and inh maps 
pia_input=original_maps(:,end);
incl_idx=1;




%% %% Get 16x16 maps for ex and in 
ex_map = reshape(original_maps(:,1:256)',16,16,length(nan_vector));
in_map = reshape(original_maps(:,257:512)',16,16,length(nan_vector));
% Get 16x16 maps for ex and in RAW
ex_map_raw = reshape(original_maps_raw(:,1:256)',16,16,length(nan_vector));
in_map_raw = reshape(original_maps_raw(:,257:512)',16,16,length(nan_vector));
%Calculate simple difference between maps
diff_map=ex_map-in_map;
%% Calculate fraction
[frac_exh abs_exh frac_inh abs_inh frac_exv abs_exv frac_inv abs_inv L23h L4h L5h pialD layer_assign] = iviv_profiles(nan_vector(incl_idx:end),str);
frac_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) frac_exv];
abs_exv_m=[zeros(length(nan_vector(incl_idx:end)),2) abs_exv];
layer_assign=[zeros(length(nan_vector(incl_idx:end)),2) layer_assign'];
frac_v=[frac_exv_m frac_inv];
frac_h=[frac_exh frac_inh];
%% Add fraction to structure 
for i=1:length(nan_vector)
    str_m(nan_vector(i)).frac_vert=frac_v(i,:);
    str_m(nan_vector(i)).frac_hori=frac_h(i,:);
end
%% Add centroid measurements
for i=1:length(nan_vector)
somax(i)=str(nan_vector(i)).somaCenter(1);
somay(i)=552-str(nan_vector(i)).somaCenter(2);
end
%IN ang centroid
[out_ang_in] = centroid_map(in_map(:,:,:),somax,pia_input,[1:cells],0);
[out_ang_inL23] = centroid_map(in_map(3:5,:,:),somax,pia_input,[1:cells],2);
[out_ang_inL4] = centroid_map(in_map(6:7,:,:),somax,pia_input,[1:cells],5);
%EX ang centroid
[out_ang_ex] = centroid_map(ex_map(:,:,:),somax,pia_input,[1:cells],0);
[out_ang_exL23] = centroid_map(ex_map(3:5,:,:),somax,pia_input,[1:cells],2);
[out_ang_exL4] = centroid_map(ex_map(6:7,:,:),somax,pia_input,[1:cells],5);
%% Add fraction to structure 
for i=1:length(nan_vector)
    str(nan_vector(i)).ang_exL23=out_ang_exL23(i,:);
    str(nan_vector(i)).ang_inL23=out_ang_inL23(i,:);
     str(nan_vector(i)).ang_exL4=out_ang_exL4(i,:);
    str(nan_vector(i)).ang_inL4=out_ang_inL4(i,:);
end
%% Add PCs

%% Align the exc and inh maps vertically
% distance between squares
grid_distance = 69;
% get the square edges
edges = 0:grid_distance:16*grid_distance;
% get the exc and inh maps
% clean_ex = reshape(cat(3,str(non_nan_cells).excMap),256,[]);
% clean_in = reshape(cat(3,str(non_nan_cells).inhMap),256,[]);
clean_ex = cell_cell(:,1:256)';
clean_in = cell_cell(:,257:512)';
clean_exraw = cell_cell_raw(:,1:256)';
clean_inraw = cell_cell_raw(:,257:512)';
pia_bins = discretize(cell_cell(:,513),edges);
% get the number of actually populated bins
bin_num = max(pia_bins);
% allocate memory for the aligned maps
aligned_maps_ex = zeros(cell_num,16+bin_num,16);
aligned_maps_in = zeros(cell_num,16+bin_num,16);
aligned_maps_exraw = zeros(cell_num,16+bin_num,16);
aligned_maps_inraw = zeros(cell_num,16+bin_num,16);

% for all the cells
for cells = 1:cell_num
    % get the cells displacement
    cell_soma = pia_bins(cells);   
    % place the map into the aligned matrix based on this displacement
    aligned_maps_ex(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_ex(:,cells),16,16);
    aligned_maps_in(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_in(:,cells),16,16);
    aligned_maps_exraw(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_exraw(:,cells),16,16);
    aligned_maps_inraw(cells,(1:16)+(1+bin_num-cell_soma),:) = reshape(clean_inraw(:,cells),16,16);
    % if the morpho is missing, skip
  
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
%% Align the exc and inh maps horizontally
clean_ex = aligned_maps_ex';
clean_in = aligned_maps_in';
% get the center of the soma
soma_centers = cat(1,str(:).somaCenter);
soma_centers = soma_centers(non_nan_cells,1)- min(soma_centers(non_nan_cells,1));
% % define the cell number
% cell_num = size(clean_ex,2);
% bin the pial vector
% pia_bins = discretize(cell_cell(:,513),bin_num);
% center_bins = discretize(soma_centers, hbin_num);
center_bins = discretize(soma_centers, edges);
% if all the cells are in the same bin, don't center
if length(unique(center_bins)) == 1
    hbin_num = 0;
else
    % get the number of actually populated bins
    hbin_num = max(center_bins);
    % allocate memory for the aligned maps
    aligned_maps_ex = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_in = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_exraw = zeros(cell_num,16+bin_num, 16+hbin_num);
    aligned_maps_inraw = zeros(cell_num,16+bin_num, 16+hbin_num);

    % for all the cells
    for cells = 1:cell_num
        % get the cells displacement
        cell_soma = center_bins(cells);   
        % place the map into the aligned matrix based on this displacement
        aligned_maps_ex(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_ex(:,cells),16+bin_num,16);
        aligned_maps_in(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_in(:,cells),16+bin_num,16);
        aligned_maps_exraw(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_exraw(:,cells),16+bin_num,16);
        aligned_maps_inraw(cells, :, (1:16)+(1+hbin_num-cell_soma)) = reshape(clean_inraw(:,cells),16+bin_num,16);
       
        % if the morpho is missing, skip
       
        
    end
end
aligned_maps_ex = reshape(aligned_maps_ex,cell_num,[]);
aligned_maps_in = reshape(aligned_maps_in,cell_num,[]);
aligned_maps_exraw = reshape(aligned_maps_exraw,cell_num,[]);
aligned_maps_inraw = reshape(aligned_maps_inraw,cell_num,[]);
%% Exctract PCs
[coeff_ex,score_ex,latent_ex,~,explained_ex,mu] = pca(aligned_maps_ex(incl_idx:end,:,:));
[coeff_in,score_in,latent_in,~,explained_in,mu] = pca(aligned_maps_in(incl_idx:end,:,:));
[coeff_exraw,score_exraw,latent_exraw,~,explained_exraw,mu] = pca(aligned_maps_exraw(incl_idx:end,:,:));
[coeff_inraw,score_inraw,latent_inraw,~,explained_inraw,mu] = pca(aligned_maps_inraw(incl_idx:end,:,:));
pc_scores=[score_ex(:,1:3) score_in(:,1:3)];
%% Add PCs to strcuture

for i=1:length(nan_vector)
    str(nan_vector(i)).PCs=pc_scores(i,:);
end
%% Clustering
[idx_input, clustering_input, leafOrder] = hca([score_ex(:,1:3) score_in(:,1:3)],0,'ward',4,pia_input,0,0.6);%call function for clustering
%% add cluster idx
for i=1:length(nan_vector)
    str(nan_vector(i)).Cluster_id=idx_input(i);
end
