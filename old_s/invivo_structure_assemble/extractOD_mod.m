function [OD] = extractOD_mod(peaks,Fit)

  %read out data/paramters from peaks & Fit structure generated via the Analysis
  %pipeline of Tobias Rose @ I:\Simon Weiler\AnalyzedData
       for k=1:size(peaks,2)
           contra_resp(k)=peaks(k).anova_p_mean_trial_contra<=0.05;
           ipsi_resp(k)=peaks(k).anova_p_mean_trial_ipsi<=0.05;
           oresp(k)=peaks(k).responder_mean_trial_anov;
           contra_tun(k)=peaks(k).Tune_Anova_maxAmpDelAve_contra;
           ipsi_tun(k)=peaks(k).Tune_Anova_maxAmpDelAve_ipsi;
           
            gOSI_ipsi(k)=1-peaks(k).ori_circvar_ipsi;
            gOSI_contra(k)=1-peaks(k).ori_circvar_contra;
            gDSI_ipsi(k)=1-peaks(k).dir_circvar_ipsi;
            gDSI_contra(k)=1-peaks(k).dir_circvar_contra;
            ODI(k)=peaks(k).ODI_final;
%            
%             mean_trial_contra(:,:,k)=peaks(k).mean_trial_contra;
%             mean_trial_ipsi(:,:,k)=peaks(k).mean_trial_ipsi;
            deltapeaks_averagetrace_contra(:,k)=peaks(k).deltapeaks_averagetrace_contra;
            deltapeaks_averagetrace_ipsi(:,k)=peaks(k).deltapeaks_averagetrace_ipsi;
%             zscore_peaks_contra(:,:,k)=peaks(k).zscore_peaks_contra;
%             zscore_peaks_ipsi(:,:,k)=peaks(k).zscore_peaks_ipsi;
%             zscore_peaks_average_contra(:,k)=peaks(k).zscore_peaks_average_contra;
%             zscore_peaks_average_ipsi(:,k)=peaks(k).zscore_peaks_average_ipsi;
%            %Reliability for each cell across all stimuli concatenated
%            cR=vertcat(peaks(k).deltatrace_trials_oris_contra{:,:});
%            iR=vertcat(peaks(k).deltatrace_trials_oris_ipsi{:,:});
%            vari_c(k)=nanmean(nanvar(cR)./nanmean(cR));
%            vari_i(k)=nanmean(nanvar(iR)./nanmean(iR));%nanmean(nanvar(cR')./nanmean(cR'))
%           
%            %Reliability of preferred ori/dir
%            idx_c= find(peaks(k).deltapeaks_averagetrace_contra==max(peaks(k).deltapeaks_averagetrace_contra));
%            idx_i= find(peaks(k).deltapeaks_averagetrace_ipsi==max(peaks(k).deltapeaks_averagetrace_ipsi));
%            vari_c_pref(k)=nanmean(nanvar(peaks(k).deltatrace_trials_oris_contra{idx_c(1),:})./nanmean(peaks(k).deltatrace_trials_oris_contra{idx_c(1),:}));
%            vari_i_pref(k)=nanmean(nanvar(peaks(k).deltatrace_trials_oris_ipsi{idx_i(1),:})./nanmean(peaks(k).deltatrace_trials_oris_ipsi{idx_i(1),:}));
%     
%            cR=[];
%            iR=[];
%            idx_c=[];
%            idx_i=[];
%            
%            
%            
%            %PSTH raw data 4s per 8 orientations
           PSTH_raw_contra{:,k}=peaks(k).deltatrace_trials_oris_contra;
           PSTH_raw_ipsi{:,k}=peaks(k).deltatrace_trials_oris_ipsi;
           
           %from Fit structure 
           prefOri_ipsi(k)=Fit(k).ipsi.PrefOri;
            prefOri_contra(k)=Fit(k).contra.PrefOri;
           prefDir_ipsi(k)=Fit(k).ipsi.PrefDir;
            prefDir_contra(k)=Fit(k).contra.PrefDir;
            oppResp_ipsi(k)=Fit(k).ipsi.OppResp;
           oppResp_contra(k)=Fit(k).contra.OppResp;
             fit_contra(:,k)=Fit(k).contra.FittedData;
             fit_ipsi(:,k)=Fit(k).ipsi.FittedData;
            sigma_ipsi(k)=Fit(k).ipsi.Sigma;
             sigma_contra(k)=Fit(k).contra.Sigma;
            error_contra(k)=Fit(k).contra.Error;
            rsquare_contra(k)=Fit(k).contra.R2;
            error_ipsi(k)=Fit(k).ipsi.Error;
            rsquare_ipsi(k)=Fit(k).ipsi.R2;
       end
      %read contra/ipsi/bino cells
      contra_only=find(contra_resp==1 & ipsi_resp==0);
      ipsi_only=find(contra_resp==0 & ipsi_resp==1);
      bino=find(contra_resp==1 & ipsi_resp==1);
      unres=find(oresp==0);
      resp=find(oresp==1);
     
      
      %contra/ipsi/bino indices
      OD.contra_only=contra_only;
      OD.ipsi_only=ipsi_only;
      OD.bino=bino;
      OD.unres=unres;
      OD.resp=resp;
      OD.oresp=oresp;
      OD.contra=contra_resp;
      OD.ipsi=ipsi_resp;
      %OD.Fit=Fit;
    
%       %OSI/DSI... based on ipsi/contra bino
     OD.gOSI=[gOSI_contra;gOSI_ipsi]' ;
     OD.gDSI=[gDSI_contra;gDSI_ipsi]' ;
     OD.ODI=ODI;
     OD.deltapeaks_averagetrace=[deltapeaks_averagetrace_contra;deltapeaks_averagetrace_ipsi]';
     OD.contra_tun=contra_tun;
     OD.ipsi_tun=ipsi_tun;
%       OD.gOSI_b=max([gOSI_ipsi(bino);gOSI_contra(bino)]);
%       OD.gDSI_c=gDSI_contra(contra_only);
%       OD.gDSI_i=gDSI_ipsi(ipsi_only);
%       OD.gDSI_b=max([gDSI_ipsi(bino);gDSI_contra(bino)]);
%       OD.ODI_c=ODI(contra_only);
%       OD.ODI_i=ODI(ipsi_only);
%       OD.ODI_b=ODI(bino);
%       OD.ODI_r=ODI(resp);
%      
%       OD.trial_contra=mean_trial_contra(:,:,contra_only);
%       OD.trial_ipsi=mean_trial_ipsi(:,:,ipsi_only);
%       OD.ave_contra=deltapeaks_averagetrace_contra(:,contra_only);
%       OD.ave_ipsi=deltapeaks_averagetrace_ipsi(:,ipsi_only);
%       OD.z_trial_contra=zscore_peaks_contra(:,:,contra_only);
%       OD.z_trial_ipsi=zscore_peaks_ipsi(:,:,ipsi_only);
%       OD.z_ave_contra=zscore_peaks_average_contra(:,contra_only);
%       OD.z_ave_ipsi=zscore_peaks_average_ipsi(:,ipsi_only);
%       %trials 
%       OD.trial_bino_c=mean_trial_contra(:,:,bino);
%       OD.trial_bino_i=mean_trial_ipsi(:,:,bino);
%       OD.ave_bino_c=deltapeaks_averagetrace_contra(:,bino);
%       OD.ave_bino_i=deltapeaks_averagetrace_ipsi(:,bino);
%       OD.z_trial_bino_c=zscore_peaks_contra(:,:,bino);
%       OD.z_trial_bino_i=zscore_peaks_ipsi(:,:,bino);
%       OD.z_ave_bino_c=zscore_peaks_average_contra(:,bino);
%       OD.z_ave_bino_i=zscore_peaks_average_ipsi(:,bino);
%       
%       %PSTH raw data 4s per 8 orientations
%       OD.PSTH_raw_contra=PSTH_raw_contra;
%       OD.PSTH_raw_ipsi=PSTH_raw_ipsi;
      
      %Fit
       OD.prefOri=[prefOri_contra;prefOri_ipsi]';
       OD.prefDir=[prefDir_contra;prefDir_ipsi]';
       OD.oppResp=[oppResp_contra;oppResp_ipsi]';
       OD.sigma=[sigma_contra;sigma_ipsi]';
       OD.fit_tuning=[fit_contra;fit_ipsi]'; 
       OD.fit_error=[error_contra;error_ipsi]';
       OD.fit_r2=[rsquare_contra; rsquare_ipsi]';
     
    
%       OD.prefDir_contra=prefDir_contra(contra_only);
%       OD.prefOri_ipsi=prefOri_ipsi(ipsi_only);
%       OD.prefOri_contra=prefOri_contra(contra_only);
%       OD.prefOri_bino_c=prefOri_contra(bino);
%       OD.prefOri_bino_i=prefOri_ipsi(bino);
%       
%       OD.prefDir_ipsi=prefDir_ipsi(ipsi_only);
%       OD.prefDir_contra=prefDir_contra(contra_only);
%       OD.prefDir_bino_c=prefDir_contra(bino);
%       OD.prefDir_bino_i=prefDir_ipsi(bino);
%       
%       OD.oppResp_ipsi=oppResp_ipsi(ipsi_only);
%       OD.oppResp_contra=oppResp_contra(contra_only);
%       OD.oppResp_bino_c=oppResp_contra(bino);
%       OD.oppResp_bino_i=oppResp_ipsi(bino);
%       
%       OD.fit_contra=fit_contra(:,contra_only);
%       OD.fit_ipsi=fit_ipsi(:,ipsi_only);
%       OD.fit_bino_c=fit_contra(:,bino);
%       OD.fit_bino_i=fit_ipsi(:,bino);
%       
%       OD.sigma_ipsi=sigma_ipsi(ipsi_only);
%       OD.sigma_contra=sigma_contra(contra_only);
%       OD.sigma_bino_c=sigma_contra(bino);
%       OD.sigma_bino_i=sigma_ipsi(bino);
%       
%       OD.error_contra=error_contra(contra_only);
%       OD.error_ipsi=error_ipsi(ipsi_only);
%       OD.error_bino_c=error_contra(bino);
%       OD.error_bino_i=error_ipsi(bino);
%       
%       OD.rsquare_contra=rsquare_contra(contra_only);
%       OD.rsquare_ipsi=rsquare_ipsi(ipsi_only);
%       OD.rsquare_bino_c=rsquare_contra(bino);
%       OD.rsquare_bino_i=rsquare_ipsi(bino);
%       
%       %Reliability 
%       OD.vari_contra=vari_c(contra_only);
%       OD.vari_ipsi=vari_i(ipsi_only);
%       OD.vari_bino_c=vari_c(bino);
%       OD.vari_bino_i=vari_i(bino);
%       %Reliability pref Ori
%       OD.varipref_contra=vari_c_pref(contra_only);
%       OD.varipref_ipsi=vari_i_pref(ipsi_only);
%       OD.varipref_bino_c=vari_c_pref(bino);
%       OD.varipref_bino_i=vari_i_pref(bino);
      
      %empty variables
      contra_resp=[];
      ipsi_resp=[];
      resp=[];
      gOSI_ipsi=[];
      gOSI_contra=[];
      gDSI_ipsi=[];
      gDSI_contra=[];
      ODI=[];
      mean_trial_contra=[];
      mean_trial_ipsi=[];
      deltapeaks_averagetrace_contra=[];
      deltapeaks_averagetrace_ipsi=[];
      zscore_peaks_contra=[];
      zscore_peaks_ipsi=[];
      zscore_peaks_average_contra=[];
      zscore_peaks_average_ipsi=[];
      
      prefOri_ipsi=[];
      prefOri_contra=[]
      prefDir_ipsi=[];
      prefDir_contra=[];
      oppResp_ipsi=[];
      oppResp_contra=[];
      fit_contra=[];
      fit_ipsi=[];
      rsquare_ipsi=[];
      rsquare_contra=[];
        error_ipsi=[];
      error_contra=[];
      contra_tun=[];
      ipsi_tun=[];
      
      PSTH_raw_ipsi={};
      PSTH_raw_contra={};
end