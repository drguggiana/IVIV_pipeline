function [od_out sftf_out spon_out] = concat_invivo(L23_PC)
for i=1:length(L23_PC)
    pia{i,:}=L23_PC(i).pial_depth;
     con_all{i,:}=L23_PC(i).OD.contra;
     ips_all{i,:}=L23_PC(i).OD.ipsi;
     ODI_all{i,:}=L23_PC(i).OD.ODI;
     prefOri_all{i,:}=L23_PC(i).OD.prefOri(:,:);
     prefDir_all{i,:}=L23_PC(i).OD.prefDir(:,:);
     OSIall{i,:}=L23_PC(i).OD.gOSI(:,:); 
     DSIall{i,:}=L23_PC(i).OD.gDSI(:,:);   
     sf_all{i,:}=L23_PC(i).SFTF.sf;
     tf_all{i,:}=L23_PC(i).SFTF.tf;
     osi_all_sftf{i,:}=L23_PC(i).SFTF.OSI;
     dsi_all_sftf{i,:}=L23_PC(i).SFTF.DSI;
     oripref_all_sftf{i,:}=L23_PC(i).SFTF.oripref;
     dirpref_all_sftf{i,:}=L23_PC(i).SFTF.dirpref;
     sad_all{i,:}=L23_PC(i).spon.sad;
     pci_all{i,:}=L23_PC(i).spon.pci;
     sftf_resp_all{i,:}=L23_PC(i).SFTF.ov_resp;
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
  OSI_all=vertcat(OSIall{:});
  DSI_all=vertcat(DSIall{:});
  %Pial depth
  pial_all=horzcat(pia{:});
  
  %Spon protocol
  sad_a=horzcat(sad_all{:});
  pci_a=vertcat(pci_all{:});
  %SFTF protocol
  sf_a=horzcat(sf_all{:});
  tf_a=horzcat(tf_all{:});
  sftfs_res_a=horzcat(sftf_resp_all{:});
  sftf_ori_a=horzcat(sftf_ori{:});
  osi_sftf=horzcat(osi_all_sftf{:});
  dsi_sftf=horzcat(dsi_all_sftf{:});
  ori_sftf=horzcat(oripref_all_sftf{:});
  dir_sftf=horzcat(dirpref_all_sftf{:});
  
  %% OD protocol
 
  for i=1:length(odi)
      if con(i)==1 & ips(i)==0;
          
          pORI_a(i)=pOSI_all(i,1);
          pDIR_a(i)=pDSI_all(i,1);
          pOSI_a(i)=OSI_all(i,1);
          pDSI_a(i)=DSI_all(i,1);
          
      elseif con(i)==0 & ips(i)==1;
          pORI_a(i)=pOSI_all(i,2);
          pDIR_a(i)=pDSI_all(i,2);
          pOSI_a(i)=OSI_all(i,2);
          pDSI_a(i)=DSI_all(i,2);
      elseif con(i)==1 & ips(i)==1;
          if odi(i)>0
              pORI_a(i)=pOSI_all(i,1);
              pDIR_a(i)=pDSI_all(i,1);
              pOSI_a(i)=OSI_all(i,1);
              pDSI_a(i)=DSI_all(i,1);
          else odi(i)<0
              pORI_a(i)=pOSI_all(i,2);
              pDIR_a(i)=pDSI_all(i,2);
              pDSI_a(i)=DSI_all(i,2);
              pOSI_a(i)=OSI_all(i,2);
          end
      else con(i)==0 & ips(i)==0
       pORI_a(i)=NaN;
       pDIR_a(i)=NaN;
       pOSI_a(i)=NaN;
       pDSI_a(i)=NaN;
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
      else 
       sf_com(i)=NaN;
       tf_com(i)=NaN;
      end
  end
  %% 
  od_out=[pOSI_a;pDSI_a ;pORI_a ;pDIR_a; ODI_res; con; ips; pial_all]';
  spon_out=[sad_a;pci_a']';
  sftf_out=[sf_com;tf_com]';
  
  
end