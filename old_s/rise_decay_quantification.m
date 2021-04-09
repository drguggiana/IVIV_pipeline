function output=rise_decay_quantification(data)
%VS 25/01/21
%VS 26/03/21 - comments, additional outputs

%input: vector of horizontal input profile
%output: 
%output.LHS_rise - spatial length between points of 10% and 90% of peak onleft hand side(LHS)
%output.RHS_rise - spatial length between points of 10% and 90% of peak on right hand side (RHS)
%output.LHS_rate - rate of rise in precent per spatial unit on LHS
%output.RHS_rate - rate of rise in precent per spatial unit on RHS
%output.disparitiyIndex - this is zero for equal rates and if one rate is far larger than the other approaches 1 (LHS > RHS) or -1 (RHS > LHS)

perc_low=0.1;
perc_high=0.9;

[data_max, index_max]=max(data);
%LHS
tempIndex=find(data(1:index_max)<(perc_low*data_max)); %LHS low
sort(tempIndex);
LHS_low=tempIndex(end)+...
    (perc_low*data_max-data(tempIndex(end)))/...
    (data(tempIndex(end)+1)-data(tempIndex(end))); %linear interpolation
tempIndex=find(data(1:index_max)<(perc_high*data_max)); %LHS high
sort(tempIndex);
LHS_high=tempIndex(end)+...
    (perc_high*data_max-data(tempIndex(end)))/...
    (data(tempIndex(end)+1)-data(tempIndex(end))); %linear interpolation
output.LHS_rise=LHS_high-LHS_low;

%RHS
tempIndex=find(data(index_max:end)<(perc_low*data_max)); %RHS low
tempIndex=index_max+sort(tempIndex)-1;
RHS_low=tempIndex(1)-...
    (perc_low*data_max-data(tempIndex(1)))/...
    (data(tempIndex(1)-1)-data(tempIndex(1))); %linear interpolation
tempIndex=find(data(index_max:end)<(perc_high*data_max)); %RHS high
tempIndex=index_max+sort(tempIndex)-1;
RHS_high=tempIndex(1)-...
    (perc_high*data_max-data(tempIndex(1)))/...
    (data(tempIndex(1)-1)-data(tempIndex(1))); %linear interpolation
output.RHS_rise=abs(RHS_high-RHS_low);

output.LHS_rate=80/output.LHS_rise;
output.RHS_rate=80/output.RHS_rise;
output.disparitiyIndex=2*output.LHS_rate/(output.LHS_rate+output.RHS_rate)-1;
end
