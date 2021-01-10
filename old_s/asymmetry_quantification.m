function output=asymmetry_quantification(data)
%VS 22/12/20
%see https://dsp.stackexchange.com/questions/15916/how-to-quantify-asymmetry-of-a-signal https://dsp.stackexchange.com/questions/15916/how-to-quantify-asymmetry-of-a-signal 
%input: vector (e.g. horizontal input profile)
%output: see below 
fsym=0.5*(data(1:end)+data(end:-1:1));
fasym=0.5*(data(1:end)-data(end:-1:1));
output.symmetricComponent=sqrt(sum(fsym.^2)); %symmetric component
output.asymmetricComponent=sqrt(sum(fasym.^2)); %asymmetric component
output.relative_asymmetry=sqrt(sum(fasym.^2))/(sqrt(sum(fsym.^2))+sqrt(sum(fasym.^2))); %relative asymmetric component

%for illustration:
% figure;hold
% plot(fsym,'red')
% plot(fasym,'blue')
% plot(data,'black')

end
