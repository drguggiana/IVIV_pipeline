function [od_out sftf_out sftf_sel sftf_pref spon_out pia_all delta_Ca fit_Ca  Ca_peak_od] = concat_invivo(L23_PC)
for i=1:length(L23_PC)
    %pia
    pia{i,:}=L23_PC(i).pial_depth;
    %OD protocol
     con_all{i,:}=L23_PC(i).OD.contra;
     ips_all{i,:}=L23_PC(i).OD.ipsi;
     ODI_all{i,:}=L23_PC(i).OD.ODI;
     prefOri_all{i,:}=L23_PC(i).OD.prefOri(:,:);
     prefDir_all{i,:}=L23_PC(i).OD.prefDir(:,:);
     OSIall{i,:}=L23_PC(i).OD.gOSI(:,:); 
     DSIall{i,:}=L23_PC(i).OD.gDSI(:,:); 
     Ca_od{i,:}=L23_PC(i).OD.deltapeaks_averagetrace(:,:);
     fit_od{i,:}=L23_PC(i).OD.fit_tuning(:,:);
     sigma{i,:}=L23_PC(i).OD.sigma(:,:);
      %SFTF protocol
     sf_all{i,:}=L23_PC(i).SFTF.sf;
     tf_all{i,:}=L23_PC(i).SFTF.tf;
     osi_all_sftf{i,:}=L23_PC(i).SFTF.OSI;
     dsi_all_sftf{i,:}=L23_PC(i).SFTF.DSI;
     oripref_all_sftf{i,:}=L23_PC(i).SFTF.oripref;
     dirpref_all_sftf{i,:}=L23_PC(i).SFTF.dirpref;
     Ca_peak_sftf{i,:}=L23_PC(i).SFTF.peakresp;
     sftf_resp_all{i,:}=L23_PC(i).SFTF.ov_resp;
     
     
     osi_sftf1{i,:}=L23_PC(i).SFTF.OSI1;
     osi_sftf2{i,:}=L23_PC(i).SFTF.OSI2;
     osi_sftf3{i,:}=L23_PC(i).SFTF.OSI3;
     osi_sftf4{i,:}=L23_PC(i).SFTF.OSI4;
     osi_sftf5{i,:}=L23_PC(i).SFTF.OSI5;
     osi_sftf6{i,:}=L23_PC(i).SFTF.OSI6;
     osi_sftf7{i,:}=L23_PC(i).SFTF.OSI7;
     osi_sftf8{i,:}=L23_PC(i).SFTF.OSI8;
     osi_sftf9{i,:}=L23_PC(i).SFTF.OSI9;
     
     dsi_sftf1{i,:}=L23_PC(i).SFTF.DSI1;
     dsi_sftf2{i,:}=L23_PC(i).SFTF.DSI2;
     dsi_sftf3{i,:}=L23_PC(i).SFTF.DSI3;
     dsi_sftf4{i,:}=L23_PC(i).SFTF.DSI4;
     dsi_sftf5{i,:}=L23_PC(i).SFTF.DSI5;
     dsi_sftf6{i,:}=L23_PC(i).SFTF.DSI6;
     dsi_sftf7{i,:}=L23_PC(i).SFTF.DSI7;
     dsi_sftf8{i,:}=L23_PC(i).SFTF.DSI8;
     dsi_sftf9{i,:}=L23_PC(i).SFTF.DSI9;
     
     ori_sftf1{i,:}=L23_PC(i).SFTF.oripref1;
     ori_sftf2{i,:}=L23_PC(i).SFTF.oripref2;
     ori_sftf3{i,:}=L23_PC(i).SFTF.oripref3;
     ori_sftf4{i,:}=L23_PC(i).SFTF.oripref4;
     ori_sftf5{i,:}=L23_PC(i).SFTF.oripref5;
     ori_sftf6{i,:}=L23_PC(i).SFTF.oripref6;
     ori_sftf7{i,:}=L23_PC(i).SFTF.oripref7;
     ori_sftf8{i,:}=L23_PC(i).SFTF.oripref8;
     ori_sftf9{i,:}=L23_PC(i).SFTF.oripref9;
     
     dir_sftf1{i,:}=L23_PC(i).SFTF.dirpref1;
      dir_sftf2{i,:}=L23_PC(i).SFTF.dirpref2;
      dir_sftf3{i,:}=L23_PC(i).SFTF.dirpref3;
      dir_sftf4{i,:}=L23_PC(i).SFTF.dirpref4;
      dir_sftf5{i,:}=L23_PC(i).SFTF.dirpref5;
      dir_sftf6{i,:}=L23_PC(i).SFTF.dirpref6;
      dir_sftf7{i,:}=L23_PC(i).SFTF.dirpref7;
      dir_sftf8{i,:}=L23_PC(i).SFTF.dirpref8;
      dir_sftf9{i,:}=L23_PC(i).SFTF.dirpref9;
      
    peak_sftf1{i,:}=L23_PC(i).SFTF.peakresp1;
    peak_sftf2{i,:}=L23_PC(i).SFTF.peakresp2;
    peak_sftf3{i,:}=L23_PC(i).SFTF.peakresp3;
    peak_sftf4{i,:}=L23_PC(i).SFTF.peakresp4;
    peak_sftf5{i,:}=L23_PC(i).SFTF.peakresp5;
    peak_sftf6{i,:}=L23_PC(i).SFTF.peakresp6;
    peak_sftf7{i,:}=L23_PC(i).SFTF.peakresp7;
    peak_sftf8{i,:}=L23_PC(i).SFTF.peakresp8;
    peak_sftf9{i,:}=L23_PC(i).SFTF.peakresp9;
    
    resp1{i,:}=L23_PC(i).SFTF.res(1,:);
    resp2{i,:}=L23_PC(i).SFTF.res(2,:);
    resp3{i,:}=L23_PC(i).SFTF.res(3,:);
    resp4{i,:}=L23_PC(i).SFTF.res(4,:);
    resp5{i,:}=L23_PC(i).SFTF.res(5,:);
    resp6{i,:}=L23_PC(i).SFTF.res(6,:);
    resp7{i,:}=L23_PC(i).SFTF.res(7,:);
    resp8{i,:}=L23_PC(i).SFTF.res(8,:);
    resp9{i,:}=L23_PC(i).SFTF.res(9,:);
    
     %Spon protocol
     sad_all{i,:}=L23_PC(i).spon.sad;
     pci_all{i,:}=L23_PC(i).spon.pci;
     
     if isempty(L23_PC(i).SFTF.Fit)==0
     sftf_ori{i,:}=[L23_PC(i).SFTF.Fit.PrefRsp];
     else
     sftf_ori{i,:}=NaN;
     end
end
  %% concatenate all
 %OD protocol
  con=horzcat(con_all{:});
  ips=horzcat(ips_all{:});
  odi=horzcat(ODI_all{:});
  pOSI_all=vertcat(prefOri_all{:});
  pDSI_all=vertcat(prefDir_all{:});
  
  temp= pOSI_all+90
temp(temp>179.999) = temp(temp>179.999) -180
pOSI_all=[];
 pOSI_all=temp;
 
  OSI_all=vertcat(OSIall{:});
  DSI_all=vertcat(DSIall{:});
  Ca_peak_od=vertcat(Ca_od{:});
  fit_res=vertcat(fit_od{:});
   tw=vertcat(sigma{:});
  %Pial depth
  pial_all=horzcat(pia{:});
  
  %Spon protocol
  sad_a=horzcat(sad_all{:});
  pci_a=vertcat(pci_all{:});
  %SFTF protocol
  sf_a=horzcat(sf_all{:});
  tf_a=horzcat(tf_all{:});
  sftfs_res_a=horzcat(sftf_resp_all{:});
  %sftf_ori_a=horzcat(sftf_ori{:});
  
  osi_sftf=horzcat(osi_all_sftf{:});
  dsi_sftf=horzcat(dsi_all_sftf{:});
  ori_sftf=horzcat(oripref_all_sftf{:});
  dir_sftf=horzcat(dirpref_all_sftf{:});
  Ca_sftf=horzcat(Ca_peak_sftf{:});
  
  osi_s=[horzcat(osi_sftf1{:});horzcat(osi_sftf2{:});horzcat(osi_sftf3{:});horzcat(osi_sftf4{:});...
      horzcat(osi_sftf5{:});horzcat(osi_sftf6{:});horzcat(osi_sftf7{:});horzcat(osi_sftf8{:});horzcat(osi_sftf9{:})]';
   dsi_s=[horzcat(dsi_sftf1{:});horzcat(dsi_sftf2{:});horzcat(dsi_sftf3{:});horzcat(dsi_sftf4{:});...
      horzcat(dsi_sftf5{:});horzcat(dsi_sftf6{:});horzcat(dsi_sftf7{:});horzcat(dsi_sftf8{:});horzcat(dsi_sftf9{:})]';
  dir_s=[horzcat(dir_sftf1{:});horzcat(dir_sftf2{:});horzcat(dir_sftf3{:});horzcat(dir_sftf4{:});...
      horzcat(dir_sftf5{:});horzcat(dir_sftf6{:});horzcat(dir_sftf7{:});horzcat(dir_sftf8{:});horzcat(dir_sftf9{:})]';
    ori_s=[horzcat(ori_sftf1{:});horzcat(ori_sftf2{:});horzcat(ori_sftf3{:});horzcat(ori_sftf4{:});...
      horzcat(ori_sftf5{:});horzcat(ori_sftf6{:});horzcat(ori_sftf7{:});horzcat(ori_sftf8{:});horzcat(ori_sftf9{:})]';
  Ca_sftf_all=[horzcat(peak_sftf1{:});horzcat(peak_sftf2{:});horzcat(peak_sftf3{:});horzcat(peak_sftf4{:});...
      horzcat(peak_sftf5{:});horzcat(peak_sftf6{:});horzcat(peak_sftf7{:});horzcat(peak_sftf8{:});horzcat(peak_sftf9{:})]';
   res_sftf_all=[horzcat(resp1{:});horzcat(resp2{:});horzcat(resp3{:});horzcat(resp4{:});...
      horzcat(resp5{:});horzcat(resp6{:});horzcat(resp7{:});horzcat(resp8{:});horzcat(resp9{:})]';
  %% OD protocol
 
  for i=1:length(odi)
      if con(i)==1 & ips(i)==0;
          
          pORI_a(i)=pOSI_all(i,1);
          pDIR_a(i)=pDSI_all(i,1);
          pOSI_a(i)=OSI_all(i,1);
          pDSI_a(i)=DSI_all(i,1);
          pCa_a(i)=max(Ca_peak_od(i,1:8));
          delta_Ca(i,:)=Ca_peak_od(i,1:8);
          fit_Ca(i,:)=fit_res(i,1:360);
            pTW(i)=tw(i,1);
      elseif con(i)==0 & ips(i)==1;
          pORI_a(i)=pOSI_all(i,2);
          pDIR_a(i)=pDSI_all(i,2);
          pOSI_a(i)=OSI_all(i,2);
          pDSI_a(i)=DSI_all(i,2);
           pCa_a(i)=max(Ca_peak_od(i,9:16));
           delta_Ca(i,:)=Ca_peak_od(i,9:16);
           fit_Ca(i,:)=fit_res(i,361:end);
           pTW(i)=tw(i,2);
           
      elseif con(i)==1 & ips(i)==1;
          if odi(i)>0
              pORI_a(i)=pOSI_all(i,1);
              pDIR_a(i)=pDSI_all(i,1);
              pOSI_a(i)=OSI_all(i,1);
              pDSI_a(i)=DSI_all(i,1);
               pCa_a(i)=max(Ca_peak_od(i,1:8));
               delta_Ca(i,:)=Ca_peak_od(i,1:8);
               fit_Ca(i,:)=fit_res(i,1:360);
                pTW(i)=tw(i,1);
          else odi(i)<0
              pORI_a(i)=pOSI_all(i,2);
              pDIR_a(i)=pDSI_all(i,2);
              pDSI_a(i)=DSI_all(i,2);
              pOSI_a(i)=OSI_all(i,2);
              pCa_a(i)=max(Ca_peak_od(i,9:16));
              delta_Ca(i,:)=Ca_peak_od(i,9:16);
              fit_Ca(i,:)=fit_res(i,361:end);
              pTW(i)=tw(i,2);
          end
      else con(i)==0 & ips(i)==0
       pORI_a(i)=NaN;
       pDIR_a(i)=NaN;
       pOSI_a(i)=NaN;
       pDSI_a(i)=NaN;
       pCa_a(i)=NaN;
       delta_Ca(i,:)=ones(8,1)*NaN;
       fit_Ca(i,:)=ones(360,1)*NaN;
         pTW(i)=NaN;
      end
  end
  %% 
  
  for i=1:length(odi)
      if con(i)==1 & ips(i)==0;
      ODI_res(i)=odi(i);
      elseif con(i)==0 & ips(i)==1;
      ODI_res(i)=odi(i);
       elseif con(i)==1 & ips(i)==1;
      ODI_res(i)=odi(i);
      else con(i)==0 & ips(i)==0
          ODI_res(i)=NaN;
      end
  end
  
      %% SFTF protocol
 
  for i=1:length(sftfs_res_a)
      if sftfs_res_a(i)==1;
          sf_com(i)=sf_a(i);
          tf_com(i)=tf_a(i);
          
          osi_pref_sftf(i)=osi_sftf(i);
          dsi_pref_sftf(i)=dsi_sftf(i);
          ori_pref_sftf(i)=ori_sftf(i);
          dir_pref_sftf(i)=dir_sftf(i);
          Ca_p_sftf(i)=Ca_sftf(i);
          
          osi_sc(i,:)=osi_s(i,:);
          dsi_sc(i,:)=dsi_s(i,:);
          ori_sc(i,:)=ori_s(i,:);
          dir_sc(i,:)=dir_s(i,:);
          Ca_sc(i,:)= Ca_sftf_all(i,:);
      else 
       sf_com(i)=NaN;
       tf_com(i)=NaN;
       
        osi_pref_sftf(i)=NaN;
          dsi_pref_sftf(i)=NaN;
          ori_pref_sftf(i)=NaN;
          dir_pref_sftf(i)=NaN;
          Ca_p_sftf(i)=NaN;
          
         osi_sc(i,:)=ones(1,9)*NaN;
          dsi_sc(i,:)=ones(1,9)*NaN;
          ori_sc(i,:)=ones(1,9)*NaN;
          dir_sc(i,:)=ones(1,9)*NaN;
          Ca_sc(i,:)=ones(1,9)*NaN;
      end
  end
  %% 
 
  %% 
  for i=1:length(osi_pref_sftf)
      if isnan(osi_pref_sftf(i))==1;
         ori_pref_sftf_f(i)=NaN;
         dir_pref_sftf_f(i)=NaN;
         Ca_pref_sftft(i)=NaN;
      else 
     ori_pref_sftf_f(i)=ori_pref_sftf(i);
           dir_pref_sftf_f(i)=dir_pref_sftf(i);
           Ca_pref_sftft(i)=Ca_p_sftf(i);
      end
  end
  %% 
  for k=1:9;
   for i=1:length(osi_s)
       if res_sftf_all(i,k)==0;
          osi_sc(i,k)=NaN; 
          dsi_sc(i,k)=NaN; 
          
       else
           osi_sc(i,k)=osi_sc(i,k); 
          dsi_sc(i,k)=dsi_sc(i,k); 
       end
   end
  end
  %% 
  for k=1:9;
   for i=1:length(osi_s)
      if isnan(osi_sc(i,k))==1;
         ori_sftf_f(i,k)=NaN;
         dir_sftf_f(i,k)=NaN;
         Ca_sc(i,k)=NaN;
      else 
    ori_sftf_f(i,k)=ori_sc(i,k);
    dir_sftf_f(i,k)=dir_sc(i,k);
    Ca_sc(i,k)=Ca_sc(i,k);
      end
   end
  end
  %% 
  for i=1:length(con)
      if con(i)==0 & ips(i)==0; 
          resp(i)=0;
      else
          resp(i)=1;
      end
  end
  %% 
  od_out=[pOSI_a;pDSI_a ;ODI_res;pORI_a ;pDIR_a; pCa_a;  pTW; resp; con; ips]';
  spon_out=[sad_a;pci_a']';
  sftf_out=[sf_com' tf_com'  1-osi_pref_sftf' 1-dsi_pref_sftf' ori_pref_sftf_f' dir_pref_sftf_f' Ca_pref_sftft' sftfs_res_a']
  sftf_sel=[1-osi_sc 1-dsi_sc]; 
  sftf_pref=[ori_sftf_f dir_sftf_f Ca_sc];
  pia_all=pial_all;
  delta_Ca=delta_Ca;
  fit_Ca=fit_Ca;
  Ca_peak_od=Ca_peak_od;
  
end