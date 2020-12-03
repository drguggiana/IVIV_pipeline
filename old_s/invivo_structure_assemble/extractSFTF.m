function [SFTF] = extractSFTF(peaks)

  %read out data/paramters from peaks & Fit structure generated via the Analysis
  %pipeline of Tobias Rose @ I:\Simon Weiler\AnalyzedData
       for k=1:size(peaks.ANOVA_p_singleTF_only,1)
          resp(:,k)=[peaks.ANOVA_p_singleSF_TF_only(k,:,1) peaks.ANOVA_p_singleSF_TF_only(k,:,2) peaks.ANOVA_p_singleSF_TF_only(k,:,3)]<0.05
       end
       
       SFTF.res=resp;
       SFTF.ov_resp=sum(resp)>=1;
end