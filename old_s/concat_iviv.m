function [od_out_iviv spon_out_iviv sftf_out_iviv sftf_out_sel_iviv sftf_out_pref_iviv] = concat_iviv(str,nan_vector)
%% Extract Oripref, Dirpref and sigmatuning for ipsi, contra and bino cell. For bino cell I use the orientation preference for the 
%more dominant eye (based on ODI)
non_nan_idx=nan_vector(1:end);
for i=1:length(non_nan_idx);
%    idx_amp_inf(i,:)=str(non_nan_idx(i)).hori_peak_pl;
if ~isempty(str(non_nan_idx(i)).Ori)==1 & isempty(str(non_nan_idx(i)).unres)==1
    %ori_both(:,i)=str(non_nan_idx(i)).ori_a;
    iv_ODI(i)=str(non_nan_idx(i)).ODI;
    if str(non_nan_idx(i)).contra==1
     oripref(i)=str(non_nan_idx(i)).Ori(1); 
     dirpref(i)=str(non_nan_idx(i)).Dir(1);
     sigma(i)=str(non_nan_idx(i)).sigma(1);
     fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(1:360);
     iv_OSI(i)=str(non_nan_idx(i)).OSI(1);
     iv_DSI(i)=str(non_nan_idx(i)).DSI(1);
     iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(1);
     iv_Ca_peak(i)=str(non_nan_idx(i)).Ca_peak(1);
    elseif str(non_nan_idx(i)).ipsi==1;
      oripref(i)=str(non_nan_idx(i)).Ori(2);
      dirpref(i)=str(non_nan_idx(i)).Dir(2);
      sigma(i)=str(non_nan_idx(i)).sigma(2);
      fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(361:end);
      iv_OSI(i)=str(non_nan_idx(i)).OSI(2);
      iv_DSI(i)=str(non_nan_idx(i)).DSI(2);
      iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(2);
      iv_Ca_peak(i)=str(non_nan_idx(i)).Ca_peak(2);
    elseif str(non_nan_idx(i)).bino==1;
         if str(non_nan_idx(i)).ODI>=0;
        oripref(i)=str(non_nan_idx(i)).Ori(1);
        dirpref(i)=str(non_nan_idx(i)).Dir(1);
        sigma(i)=str(non_nan_idx(i)).sigma(1);
        fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(1:360);
        iv_OSI(i)=str(non_nan_idx(i)).OSI(1);
        iv_DSI(i)=str(non_nan_idx(i)).DSI(1);
        iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(1);
        iv_Ca_peak(i)=str(non_nan_idx(i)).Ca_peak(1);
         else str(non_nan_idx(i)).ODI<0;
            oripref(i)=str(non_nan_idx(i)).Ori(2);
        dirpref(i)=str(non_nan_idx(i)).Dir(2);
        sigma(i)=str(non_nan_idx(i)).sigma(2);
        fit_resp(i,:)=str(non_nan_idx(i)).fit_resp(361:end);
        iv_OSI(i)=str(non_nan_idx(i)).OSI(2);
        iv_DSI(i)=str(non_nan_idx(i)).DSI(2);
        iv_Ca(i)=str(non_nan_idx(i)).Ca_sum(2);
        iv_Ca_peak(i)=str(non_nan_idx(i)).Ca_peak(2);
    end
 else str(non_nan_idx(i)).unres==1;
        oripref(i)=NaN;
        dirpref(i)=NaN;
        sigma(i)=NaN;
        fit_resp(i,:)=ones(360,1)*NaN;
        iv_ODI(i)=NaN;
        iv_OSI(i)=NaN;
        iv_DSI(i)=NaN;
        iv_Ca(i)=NaN;
        iv_Ca_peak(i)=NaN;
        %ori_both(:,i)=ones(1,2)*NaN;
 end

else
oripref(i)=NaN;
dirpref(i)=NaN;
sigma(i)=NaN;
fit_resp(i,:)=ones(360,1)*NaN;
iv_ODI(i)=NaN;
iv_OSI(i)=NaN;
iv_DSI(i)=NaN;
iv_Ca(i)=NaN;
iv_Ca_peak(i)=NaN;
end
end
%% Extract other in vivo paramters: PCI and spontenaeous events 
for i=1:length(non_nan_idx);   
if ~isempty(str(non_nan_idx(i)).pci)==1
iv_spon(i)=str(non_nan_idx(i)).sad;
iv_popcop(i)=str(non_nan_idx(i)).pci;
iv_pia_input(i)=str(non_nan_idx(i)).pia_invivo;
else
iv_spon(i)=NaN;    
iv_popcop(i)=NaN; 
iv_pia_input(i)=NaN;
end
end
% %% %OD decompistion
% for i=1:length(non_nan_idx)
%     if ~isempty(str(non_nan_idx(i)).iv_OD_decom)==1
%         od_decomp(i,:)=str(non_nan_idx(i)).iv_OD_decom;
%     else
%         od_decomp(i,:)=ones(1,24)*NaN;
%     end
% end
% % SF decomposition
% for i=1:length(non_nan_idx)
%     if ~isempty(str(non_nan_idx(i)).iv_SFTF_decom)==1
%         sftf_decomp(i,:)=str(non_nan_idx(i)).iv_SFTF_decom;
%     else
%         sftf_decomp(i,:)=ones(1,24)*NaN;
%     end
% end
%% SFTF both eyes
for i=1:length(non_nan_idx)
    if ~isempty(str(non_nan_idx(i)).sftf_resp)==1
        iv_SF(i)=str(non_nan_idx(i)).SF;
        iv_TF(i)=str(non_nan_idx(i)).TF;
        oripref_sftf(i)=str(non_nan_idx(i)).Ori_sftf;
        dirpref_sftf(i)=str(non_nan_idx(i)).Dir_sftf;
        Ca_peak_sftf(i)=str(non_nan_idx(i)).Ca_peak_sftf;
        oripref_sftf_all(i,:)=str(non_nan_idx(i)).Ori_sftf_all;
        dirpref_sftf_all(i,:)=str(non_nan_idx(i)).Dir_sftf_all;
        osi_sftf(i)=str(non_nan_idx(i)).OSI_sftf;
        dsi_sftf(i)=str(non_nan_idx(i)).DSI_sftf;
        osi_sftf_all(i,:)=str(non_nan_idx(i)).OSI_sftf_all;
        dsi_sftf_all(i,:)=str(non_nan_idx(i)).DSI_sftf_all;
        Ca_peak_sftf_all(i,:)=str(non_nan_idx(i)).Ca_sftf_all;
        if sum(str(non_nan_idx(i)).sftf_resp)==0;
            iv_resp2(i)=0;
        else sum(str(non_nan_idx(i)).sftf_resp)>0;
            iv_resp2(i)=1;
        end
    else
        iv_SF(i)=NaN;
        iv_TF(i)=NaN;  
        iv_resp2(i)=NaN;
        oripref_sftf(i)=NaN;
        dirpref_sftf(i)=NaN;
        osi_sftf(i)=NaN;
        dsi_sftf(i)=NaN;
        Ca_peak_sftf(i)=NaN;
        osi_sftf_all(i,:)=ones(1,9)*NaN;
        dsi_sftf_all(i,:)=ones(1,9)*NaN;
        oripref_sftf_all(i,:)=ones(1,9)*NaN;
        dirpref_sftf_all(i,:)=ones(1,9)*NaN;
        Ca_peak_sftf_all(i,:)=ones(1,9)*NaN;
    end
end
%% Assemble out
od_out_iviv = [iv_OSI' iv_DSI' iv_ODI'  oripref' dirpref' iv_Ca' sigma' iv_Ca_peak' iv_pia_input'];
spon_out_iviv = [iv_spon' iv_popcop'];
sftf_out_iviv = [iv_resp2' iv_SF' iv_TF' osi_sftf' dsi_sftf' oripref_sftf' dirpref_sftf' Ca_peak_sftf']
sftf_out_sel_iviv=[osi_sftf_all dsi_sftf_all];
sftf_out_pref_iviv=[oripref_sftf_all dirpref_sftf_all Ca_peak_sftf_all];

end