 
%% LOAD in vivo data structure L23PC data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file
%% save location
adata_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\invivo_only_str'

%% Excel sheet directory
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\L23PC_project_all_invivo.xlsx';%directory where excel batch file is located;change accordingly
%%   Read out excel info
batchopt          = parseExperimentsXls_L23(ExpXls);%calls the nested function parseExperimentsXls_dLGN and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
%% 

adder=1;%counting variable
for i=1:length(L23_PC)
   temp=[]
    temp=L23_PC(i).ivivROI;
    
   k=[];
    cellname=cell(length(temp),1);
    cellid=cell(length(temp),1);
    binoexp=cell(length(temp),1);
   sftfexp=[];
   sponexp=cell(length(temp),1);
  for k=1:length(temp)
   cellname{k,:}=[char(batchopt.exp_invitro{i}) num2str(batchopt.invitro_rec{i}(k),'%04.f')];
   cellid{k,:}=batchopt.cellID{i}(k);
    binoexp(k,:)=batchopt.binoexp_ids(i);
    sftfexp(k,:)=batchopt.sftfexp_ids{i}(1);
    sponexp(k,:)=batchopt.sponexp_ids(i);
  end 
  
  tr(adder).cellname=cellname
  tr(adder).cellID=cellid;
  tr(adder).binoexp=binoexp;
  tr(adder).sftfexp=sftfexp;
  tr(adder).sponexp=sponexp;
adder=adder+1;
end

%% get the 77 in vivo invitor cells
cellnames_all=cat(1,tr(:).cellname);
cellids_all=cell2mat(cat(1,tr(:).cellID));
binoexp_all=cell2mat(cat(1,tr(:).binoexp));
sftfexp_all=cat(1,tr(:).sftfexp);
sponexp_all=cell2mat(cat(1,tr(:).sponexp));
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
  %% Adding cellnames ids etc
  for i=1:length(cellnames_all)
  invivo_struct(i).cellnames=cellnames_all(i);
  invivo_struct(i).cellids=cellids_all(i);
  invivo_struct(i).ODexpID=binoexp_all(i);
  invivo_struct(i).SFTFexpID=sftfexp_all(i);
  invivo_struct(i).iviv=1;
  end
  %% Adding ipsi, contra, bino, unres
  for i=1:length(cellnames_all)
    if contra(i)==1 & ipsi(i)==0
        invivo_struct(i).contra=1;
    elseif contra(i)==0 & ipsi(i)==1
   invivo_struct(i).ipsi=1;
    elseif contra(i)==1 & ipsi(i)==1
     invivo_struct(i).bino=1;   
    else contra(i)==0 & ipsi(i)==0
        invivo_struct(i).unres=1;       
    end
end
  %% ADD invivo features to strcutre, for OD protocol: first column contra, second column ipsi
   str_m=invivo_struct;
aiviv=1:length(cellnames_all);
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
aiviv=1:length(cellnames_all);
for i=1:length(aiviv)
    str_m(aiviv(i)).c_tun=contra_tun(i);
    str_m(aiviv(i)).i_tun=ipsi_tun(i);
    str_m(aiviv(i)).error_fit=error_fit(i,:);
    str_m(aiviv(i)).r2_fit=r2_fit(i,:);
    
end
%% 
invivo_struct=str_m;
%% 
[od_out_iviv spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(invivo_struct,1:length(cellnames_all));

%% 
nan_vector=1:length(cellnames_all);
for i=1:length(nan_vector)
   invivo_struct(nan_vector(i)).OSIpref=od_out_iviv(i,1);
   invivo_struct(nan_vector(i)).DSIpref=od_out_iviv(i,2);
    invivo_struct(nan_vector(i)).ODIpref=od_out_iviv(i,3);
    invivo_struct(nan_vector(i)).ORIpref=od_out_iviv(i,4);
    invivo_struct(nan_vector(i)).DIRpref=od_out_iviv(i,5);
    invivo_struct(nan_vector(i)).Casumpref=od_out_iviv(i,6);
    invivo_struct(nan_vector(i)).Sigmapref=od_out_iviv(i,7);
    invivo_struct(nan_vector(i)).Capeakpref=od_out_iviv(i,8);
    invivo_struct(nan_vector(i)).error_pref=od_out_iviv(i,10);
    invivo_struct(nan_vector(i)).r2_pref=od_out_iviv(i,11);
    invivo_struct(nan_vector(i)).tun_pref=od_out_iviv(i,12);
end
%% Save
  cd(adata_dir);
        FileName=['invivo_struct'];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName,'invivo_struct','-v7.3');
        disp('FILE SAVED');
